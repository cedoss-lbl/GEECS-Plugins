# -*- coding: utf-8 -*-
"""
Created on Fri Feb  3 10:44:15 2023

@author: ReiniervanMourik
"""

#%% init
from __future__ import annotations

from math import sqrt, cos, sin, pi
from typing import Optional, TYPE_CHECKING, Annotated

from itertools import product

import warnings

import numpy as np

from scipy.sparse import csr_array
from scipy.sparse.linalg import lsqr
from scipy.interpolate import RegularGridInterpolator

from skimage.restoration import unwrap_phase

from pint import UnitRegistry
if TYPE_CHECKING:
    from pint import Quantity
    SpatialFrequencyQuantity = Annotated[Quantity, '[length]**-1]']
    LengthQuantity = Annotated[Quantity, '[length]']

ureg = UnitRegistry()
Q_ = ureg.Quantity

from collections import namedtuple
SpatialFrequencyCoordinates = namedtuple('SpatialFrequencyCoordinates', ['nu_x', 'nu_y'])

#%% PhasicsImageAnalyzer class

class PhasicsImageAnalyzer:
    """ An engine that can analyze a Phasics image and return its phase map.

    General usage, with `img` an image from the HTU Gasjet Phasics camera 
       pia = PhasicsImageAnalyzer()
       phase_map = pia.calculate_phase_map(pia.crop_image(img))

              
    Methods
    -------
    calculate_phase_map(img)
        runs all steps of the phase map reconstruction from a Phasics image img.

    """

    CAMERA_RESOLUTION = Q_(4.74, 'micrometer')
    GRATING_CAMERA_DISTANCE = Q_(0.841, 'millimeter')
    GRATING_PITCH = Q_(19.0505, 'micrometer')
    CAMERA_TILT = Q_(0.529878, 'radians')
    
    @property
    def diffraction_spot_centers(self):
        """ Returns the 4 diffraction spot centers in the Fourier transform
        
        Besides the spot at the origin of the FT (which represents slow 
        variations across the image), there are 8 1st order diffraction spots 
        arranged in a square around the origin, tilted by CAMERA_TILT with 
        half-width 1 / GRATING_PITCH

        Somewhat arbitrarily, I've picked the spots to be
        0. The top right corner, which is sqrt(2) / GRATING_PITCH away at 
           CAMERA_TILT + 45deg 
        1. The midpoint of the right edge of the square, at CAMERA_TILT
        2. The lower right corner, again sqrt(2) / GRATING_PITCH away, at
           CAMERA_TILT - 45deg
        3. The midpoint of the bottom edge of the square, i.e. CAMERA_TILT - 90deg.

        """
        return [SpatialFrequencyCoordinates(distance_multiplier * 1 / self.GRATING_PITCH * cos(self.CAMERA_TILT + angle), 
                                            distance_multiplier * 1 / self.GRATING_PITCH * sin(self.CAMERA_TILT + angle)
                                           )
                for distance_multiplier, angle in [(sqrt(2), pi/4), (1.0, 0.0), (sqrt(2), -pi/4), (1.0, -pi/2)]
               ]

    def __init__(self,
                 reconstruction_method = 'baffou',
                 diffraction_spot_crop_radius: Optional[Quantity] = None
                ):
        """ 
        Parameters
        ----------
        reconstruction_method : 'baffou' or 'velghe'
            which method to use for the step of combining the various spot FTs
            into a final phase map.
            'baffou' recovers phase map gradients in different directions and then
            integrates them
            'velghe' (not yet implemented) solves for the FT of the phase map

        diffraction_spot_crop_radius : Quantity with [length]^-1 units
            radius of disc around spot center to use for each spot's FT.
            If None, find the maximum radius that causes no overlap.

        """

        self.reconstruction_method = reconstruction_method
        self.diffraction_spot_crop_radius = diffraction_spot_crop_radius

    
    def _fourier_transform(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ Takes the fourier transform of an image and shifts it.

        Returns
        -------
        IMG : np.ndarray
            Fourier-transformed image, with freq = 0,0 in the middle
        freq_x : np.ndarray
            array of frequencies in the axis=1 direction
        freq_y : np.ndarray
            array of frequencies in the axis=0 direction
        
        """
        fftshift = ureg.wraps('=A', ('=A', None))(np.fft.fftshift)
        self.IMG = np.fft.fftshift(np.fft.fft2(self.img))
        self.freq_x = fftshift(np.fft.fftfreq(self.shape[1], d=self.CAMERA_RESOLUTION))
        self.freq_y = fftshift(np.fft.fftfreq(self.shape[0], d=self.CAMERA_RESOLUTION))
        
        return self.IMG, self.freq_x, self.freq_y
    
    
    def locate_diffraction_spot(self, 
                                center0: SpatialFrequencyCoordinates, 
                                search_radius: SpatialFrequencyQuantity = Q_(1.2, 'mm^-1')
                               ) -> SpatialFrequencyCoordinates:
        """ Utility function to refine location of diffraction spot center

        Uses 2D parabolic interpolation to find peak in the Fourier transform 
        of the image saved in this instance, near a given center0
        
        Parameters
        ----------
        center0 : tuple[SpatialFrequencyQuantity, SpatialFrequencyQuantity]
            x and y coordinates (not row, column) of approximate center
        search_radius : SpatialFrequencyQuantity
            how far around center0 to look for the maximum.

        Returns
        -------
        center : tuple[SpatialFrequencyQuantity, SpatialFrequencyQuantity]
            x and y coordinates of refined center, in [length]^-1 units.

        """

        # first find the maximum value of IMG near center0
        FREQ_X, FREQ_Y = np.meshgrid(self.freq_x, self.freq_y)
        IMG_magn = np.abs(self.IMG)
        Z = IMG_magn * (np.square(FREQ_X - center0.nu_x) + np.square(FREQ_Y - center0.nu_y) <= search_radius**2)
        y0, x0 = np.unravel_index(np.argmax(Z.flatten()), Z.shape)

        # then use parabolic interpolation 1 pixel around x0, y0
        # We fit a 2d quadratic function to the 9 pixels around and including
        # x0, y0
        #   f(dx, dy) = coeffs . (1, dx, dy, dx^2, dx*dy, dy^2)
        # We construct matrix A where rows are (1, dx, dy, dx^2, dx*dy, dy^2) for 
        # (dx, dy) = [(0, 0), (-1, 0), (+1, 0), (0, -1), (0, +1), (-1, -1), (-1, +1), (+1, -1), (+1, +1)]
        # Furthermore, the matrix W is diag([1, 1, 1, 1, 1, 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2)]) 
        # to downweight the corners of the 3x3 square.
        # then to solve the least squares equation A . coeffs = z, with z the values
        # of IMG at (x0 + dx, y0 + dy) for the above (dx, dy), we construct
        # B = (A.T . W . A)^-1 . A.T . W
        B = [[(15 + 2*sqrt(2))/31, (8 - sqrt(2))/31, (8 - sqrt(2))/31, (8 - sqrt(2))/31, (8 - sqrt(2))/31, (-8 + sqrt(2))/62, (-8 + sqrt(2))/62, (-8 + sqrt(2))/62, (-8 + sqrt(2))/62], [0, 1/2 - 1/sqrt(2), -1/2 + 1/sqrt(2), 0, 0, (-2 + sqrt(2))/4, (-2 + sqrt(2))/4, (2 - sqrt(2))/4, (2 - sqrt(2))/4], [0, 0, 0, 1/2 - 1/sqrt(2), -1/2 + 1/sqrt(2), (-2 + sqrt(2))/4, (2 - sqrt(2))/4, (-2 + sqrt(2))/4, (2 - sqrt(2))/4], [(-7 - 3*sqrt(2))/31, (7 + 3*sqrt(2))/62, (7 + 3*sqrt(2))/62, (3*(-8 + sqrt(2)))/62, (3*(-8 + sqrt(2)))/62, (-3*(-8 + sqrt(2)))/124, (-3*(-8 + sqrt(2)))/124, (-3*(-8 + sqrt(2)))/124, (-3*(-8 + sqrt(2)))/124], [0, 0, 0, 0, 0, 1/4, -1/4, -1/4, 1/4], [(-7 - 3*sqrt(2))/31, (3*(-8 + sqrt(2)))/62, (3*(-8 + sqrt(2)))/62, (7 + 3*sqrt(2))/62, (7 + 3*sqrt(2))/62, (-3*(-8 + sqrt(2)))/124, (-3*(-8 + sqrt(2)))/124, (-3*(-8 + sqrt(2)))/124, (-3*(-8 + sqrt(2)))/124]]
        # and solve coeffs = B . z
        z = [IMG_magn[y0 + dy, x0 + dx] 
             for (dx, dy) in [(0, 0), (-1, 0), (+1, 0), (0, -1), (0, +1), (-1, -1), (-1, +1), (+1, -1), (+1, +1)]
            ]
        c_1, c_x, c_y, c_x2, c_xy, c_y2 = c = np.dot(B, z)
        # finally, we solve df/ddx = 0 and df/ddy = 0
        dx, dy = np.array([ c_xy * c_y - 2 * c_x * c_y2, 
                            c_xy * c_x - 2 * c_y * c_x2
                         ]) / (c_xy**2 - 4 * c_x2 * c_y2)

        # and we return the spatial frequencies corresponding to x0 + dx and y0 + dy
        return SpatialFrequencyCoordinates(np.interp(x0 + dx,  np.arange(len(self.freq_x)), self.freq_x),
                                           np.interp(y0 + dy,  np.arange(len(self.freq_y)), self.freq_y),
                                          ) 

    def _set_diffraction_spot_crop_radius(self) -> SpatialFrequencyQuantity:
        """ Calculate minimum distance between any two centers
        """
        self.diffraction_spot_crop_radius = np.sqrt( 
            min((center2.nu_x - center1.nu_x)**2 + (center2.nu_y - center1.nu_y)**2
                for (i1, center1), (i2, center2) in product(enumerate(self.diffraction_spot_centers), enumerate(self.diffraction_spot_centers))
                if i1 != i2
               )
        ) / 2
        
        return self.diffraction_spot_crop_radius
  


    
    def _crop_and_center_diffraction_spots(self) -> list[np.ndarray]:

        # if diffraction spot crop radius not given, infer it
        if self.diffraction_spot_crop_radius is None:
            self._set_diffraction_spot_crop_radius()
              
        def _crop_and_center_diffraction_spot(center: SpatialFrequencyCoordinates) -> np.ndarray:
            """ Crop an area of freq space around diffraction spot center, and translate it
                to the middle of the image.
            
            Parameters
            ----------
            IMG : np.ndarray
                Fourier-transformed QWLSI image.
            center : SpatialFrequencyCoordinates
                location of diffraction spot center
    
            Returns
            -------
            G_i : np.ndarray
                image of same size as IMG, with diffraction spot in center and all
                information not related to this spot cropped away.
    
            """
            
            NU_X, NU_Y = np.meshgrid(self.freq_x, self.freq_y)
            IMG_cropped = self.IMG * ((np.square(NU_X - center.nu_x) + np.square(NU_Y - center.nu_y)) < self.diffraction_spot_crop_radius**2)

            # RegularGridInterpolator will strip the units off the grid and values, 
            # so make sure they are all handled correctly.
            @ureg.wraps('=A', ('=A', '=B', '=B', '=B', '=B'))
            def recenter_IMG_cropped_ua(IMG_cropped, freq_x, freq_y, nu_x, nu_y):
                # create a interpolator of IMG_cropped on a coordinate system that has
                # this diffraction spot center as its origin.
                IMG_cropped_interpolator = RegularGridInterpolator((freq_x - nu_x, freq_y - nu_y), IMG_cropped.T, 
                                                                   method='linear', bounds_error=False, fill_value=0.0
                                                                  )
                # evaluate the interpolator at self.freq_x x self.freq_y 
                NU_X, NU_Y = np.meshgrid(freq_x, freq_y)
                return IMG_cropped_interpolator(np.stack([NU_X.flatten(), NU_Y.flatten()], axis=1)).reshape(IMG_cropped.shape)

            return recenter_IMG_cropped_ua(IMG_cropped, self.freq_x, self.freq_y, center.nu_x, center.nu_y)

        self.diffraction_spot_IMGs = [_crop_and_center_diffraction_spot(center)
                                      for center in self.diffraction_spot_centers
                                     ]
        
        return self.diffraction_spot_IMGs


    def _reconstruct_wavefront_gradients_from_cropped_centered_diffraction_FTs(self) -> list[np.ndarray]:
        """ Calculate the angle (argument) of the inverse FT of a diffraction 
            spot image.
            
            The wavefront gradient in a specific direction is in the argument 
            of the inverse FT. 
        
            In particular, these wavefront gradient maps represent
                nu . grad(W),
            where W is the wavefront, i.e. relative optical distance (in [length] 
            units)
            
            Returns
            -------
            wavefront_gradient_maps : list of Quantity np.ndarray
                has units [length]**-1

        """

        self.wavefront_gradients = [unwrap_phase(np.angle(np.fft.ifft2(np.fft.ifftshift(IMG.m))))
                                    / (2 * np.pi * self.GRATING_CAMERA_DISTANCE)
                                    for IMG in self.diffraction_spot_IMGs
                                   ] 

        return self.wavefront_gradients



    def _integrate_gradient_maps(self) -> np.ndarray:
        """ Calculates the phase map from gradients in different directions.
        
        More precisely, this finds the phase map whose gradients in different directions
        most closely match the given gradient maps. It constructs a matrix representing 
        finite differences in each gradient direction for each pixel and solves the least
        squares matrix equation.

        Returns
        -------
        wavefront : np.ndarray

        """
        
        # Construct sparse matrix specifying the gradients in each direction for each pixel
        # for the equation A.x = b. 
        # Each row of A, and its corresponding value in b, represents the gradient in one particular
        # direction for one particular pixel of the phase map. The row is constructed by taking a 2D
        # image of zeros except for a few finite difference coefficients around a specific pixel, then
        # flattening it. 

        m = 0
        data = []
        row_ind = []
        col_ind = []

        b = []

        def to_flattened_index(i, j):
            return i * self.shape[1] + j

        for center, wavefront_gradient in zip(self.diffraction_spot_centers, self.wavefront_gradients):
            # finite_difference_coefficients = np.array(
            #     [[ 0.0,               -center.nu_y / 2,       0.0       ],
            #      [ -center.nu_x / 2,         0.0       center.nu_x / 2  ],
            #      [ 0.0                 center.nu_y / 2,       0.0       ]
            #     ]
            # ) / self.CAMERA_RESOLUTION

            g_x = (center.nu_x / self.CAMERA_RESOLUTION).m_as('m^-2')
            g_y = (center.nu_y / self.CAMERA_RESOLUTION).m_as('m^-2')
            wg = wavefront_gradient.m_as('m^-1')

            for i in range(0, self.shape[0]):
                for j in range(0, self.shape[1]):
                    row_ind.extend([m, m])
                    if j == 0:
                        col_ind.extend([to_flattened_index(i, j), to_flattened_index(i, j + 1)])
                        data.extend([-g_x, g_x])
                    elif j == self.shape[1] - 1:
                        col_ind.extend([to_flattened_index(i, j - 1), to_flattened_index(i, j)])
                        data.extend([-g_x, g_x])
                    else:
                        col_ind.extend([to_flattened_index(i, j - 1), to_flattened_index(i, j + 1)])
                        data.extend([-g_x / 2 , g_x / 2])

                    row_ind.extend([m, m])
                    if i == 0:
                        col_ind.extend([to_flattened_index(i, j), to_flattened_index(i + 1, j)])
                        data.extend([-g_y, g_y])
                    elif i == self.shape[0] - 1:
                        col_ind.extend([to_flattened_index(i - 1, j), to_flattened_index(i, j)])
                        data.extend([-g_y, g_y])
                    else:
                        col_ind.extend([to_flattened_index(i - 1, j), to_flattened_index(i + 1, j)])
                        data.extend([-g_y / 2 , g_y / 2])

                    b.append(wg[i, j])
                    m += 1


        # The least squares loss is invariant under adding a constant to the entire phase map, so add a row to the list of 
        # equations requiring that the upper left pixel's phase = 0

        data.append(np.sqrt(g_x**2 + g_y**2))  # the coefficient doesn't matter, but it's good if it's in the same order of magnitude as the rest of the coefficients 
        row_ind.append(m)
        col_ind.append(to_flattened_index(0, 0))
        b.append(0.0)

        A = csr_array((data, 
                       (row_ind, col_ind)
                      ), shape=(len(self.wavefront_gradients) * self.shape[0] * self.shape[1] + 1, 
                                self.shape[0] * self.shape[1]
                               )
                      )

        # solve the linear equation in the least squares sense.
        self.wavefront = Q_(lsqr(A, b)[0].reshape(self.shape), 'm')

        return self.wavefront

    
    def _reconstruct_wavefront_FT_from_gradient_FTs(self):
        """ Fit wavefront whose FT best corresponds to its gradient FTs
        
        From the method in 
            Velghe, Sabrina, Jérôme Primot, Nicolas Guérineau, Mathieu Cohen, 
            and Benoit Wattellier. "Wave-Front Reconstruction from Multidirectional 
            Phase Derivatives Generated by Multilateral Shearing Interferometers." 
            Optics Letters 30, no. 3 (February 1, 2005): 245. 
            https://doi.org/10.1364/OL.30.000245.

        This method solves for FT(W), where W is the wavefront, by minimizing
        the error between each gradient wavefront FT G_j and the gradient of W
        in the Fourier domain, 2*pi*i*u_j*W, where u_j is the conjugate of 
        the spatial coordinate. 
        
        The wavefront is then simply the inverse FT of the best fit.

        """
        
        # Calculate the u_j. These are calculated as  nu . v, where nu is the
        # coordinates in the fourier domain, and v is the vector of the 
        # gradient in the spatial domain.
        NU_X, NU_Y = np.meshgrid(self.freq_x, self.freq_y)
        U = [  NU_X * center.nu_x + NU_Y * center.nu_y
             for center in self.diffraction_spot_centers
            ]

        # calculate FTs of wavefront gradient maps
        # G will have units of [length]^-1
        @ureg.wraps('=A', '=A')
        def fft2_with_shift_ua(wgm):
            return np.fft.fftshift(np.fft.fft2(wgm))
        G = [fft2_with_shift_ua(wgm) for wgm in self.wavefront_gradients]

        # Solve for FT(W)_e, the estimate of the FT of the wavefront.
        # suppress divide by 0 error (if centers have integer row and column, 
        # their corresponding u will have a 0.0)
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', "invalid value encountered in divide")
            
            W_ft = (-1j/(2*np.pi) * sum(u * g for u, g in zip(U, G)) 
                                  / sum(np.square(u) for u in U)
                   )

        @ureg.wraps('=A', '=A')
        def ifft2_with_shift_mean_zero_ua(W_ft):
            # W_ft has axes -Ny..Ny x -Ny..Ny. For the ifft2, shift it back to 
            # 0..2*Ny x 0..2*Ny
            W_ft = np.fft.ifftshift(W_ft)
            # after shifting, the value corresponding to freq_x = freq_y = 0 should
            # be at 0, 0, and its value should be NaN because of the division by 
            # zero in estimating W_ft. 
            assert np.isnan(W_ft[0, 0])
            # setting it to 0 ensures that the mean of the phase map is 0. 
            W_ft[0, 0] = 0.0
            
            W = np.fft.ifft2(W_ft)
            
            return W.real
        
        W = ifft2_with_shift_mean_zero_ua(W_ft)
        
        self.wavefront = W.to('nm')

        return self.wavefront
    

    
    def calculate_wavefront(self, img: np.ndarray) -> np.ndarray:
        """ Analyze cropped Phasics quadriwave shearing image
        
        Takes a cropped image and runs the full algorithm on it to obtain the
        reconstructed phase map. 

        Parameters
        ----------
        img : np.ndarray
            Phasics image, already cropped to region of interest.

        Returns
        -------
        phase map : np.ndarray
            reconstructed phase in radians.

        """


        # take 2D Fourier transform of image, shifted so that freq = 0, 0 is 
        # in the middle of the image.
        self.img = img
        self.shape = img.shape
        self._fourier_transform()

        # get cropped and centered FTs of each diffraction spot
        self._crop_and_center_diffraction_spots()

        self._reconstruct_wavefront_gradients_from_cropped_centered_diffraction_FTs()

        # reconstruct wavefront from diffraction spot FTs
        if self.reconstruction_method == 'baffou':
            W = self._integrate_gradient_maps()

        elif self.reconstruction_method == 'velghe':
            W = self._reconstruct_wavefront_FT_from_gradient_FTs()

        else:
            raise ValueError(f"Unknown reconstruction method: {self.reconstruction_method}")

        return W


    def calculate_phase_map(self, img: np.ndarray, wavelength=Q_(800, 'nm')):
        self.calculate_wavefront(img)
        phase_map = Q_(2 * np.pi, 'radian') * self.wavefront / wavelength
        return phase_map.to_base_units()

