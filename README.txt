Point Spread Function Estimation from a Random Target
======================================================
Version 1.3 - December 26, 2011 (see list of BUGS)

by    Mauricio Delbracio <mdelbra@gmail.com>
      Pablo Muse
      Andres Almansa


Introduction
-----------
Extrinsic image blur can be observed when the camera's focal distance is not
correctly adjusted by the user, when the objects in the scene  appear at
different depths, or when the relative motion between the  camera and the
scen is faster than the shutter speed (motion blur). Besides these sources
of blur, even in ideal acquisition conditions  there is a permanent intrinsic
physical camera blur due to light diffraction, sensor resolution, lens
aberration, and anti-aliasing  filters. Our goal here is to accurately
estimate the Point Spread  Function - PSF, that models the intrinsic camera
blur. This function can be locally interpreted as the response of the camera
to a point light source.

In [1] we presented a theoretical study proving that the sub-pixel PSF
estimation problem is well-posed even with a single image capture, as long as
the captured scene is well chosen. Indeed, theoretical bounds show that a
near-optimal accuracy can be achieved by taking a single snapshot of a
calibration pattern mimicking a Bernoulli(0.5) white noise. We first use an
algorithm to accurately estimate the pattern position and its illumination
conditions. This allows for accurate geometric registration and radiometric
correction; Once these procedures have been applied, the local PSF can be
directly computed by inverting a linear system that is well-posed and
consequently its inversion does not require any regularization or prior model.



 Files
 -----
 COPYING
 Makefile
 README.txt
 VERSION
 lsd.c
 lsd.h
 thin_plates.c
 thin_plates.h
 image.c
 image.h
 io_pgm.c
 io_phm.h
 detect_pattern.c
 detect_pattern.h
 psf_estim.c
 psf_estim.h
 psf_estim_main.c
 nnls.c
 nnls.h
 pattern_noise.pgm
 img_example.pgm
 doxygen.config


Requirements
------------
- The fftw3 header and libraries are required on the system for
compilation and execution. See http://www.fftw.org/

- The cblas header and libraries are required on the system for
compilation and execution.

- The lapack library is required on the system for
compilation and execution.


Compilation
-----------
Simply use the provided Makefile, with the command `make`. You need to set
the directory where the libraries: ffw3, cblas and lapack have the respective
header and libraries files.

Running
-------

Usage: ./psf_estim [options] <input file> <output file>

Only PGM 16/8 bits images are supported.

Options:
  -s <number>   The super-resolution factor, positive integer  (default 4)
  -k <number>   PSF support size (default 4s+1)
  -o <file>     Estimated PSF saved to a 8bits PGM image
  -p <file>     Noise Pattern PGM image
  -d <file>     Save the input image with detected pattern marked (24bits PPM)
  -t <0,1,2>    Least Squares solver:
                     0 - Least Squares
                     1 - Least Squares + thresholding  (default)
                     2 - Non-negative least Squares (slowest)


Parameter Explanation

-s 'number' : The superresolution factor, i.e. how many additional samples per
              observed pixel will be computed. (default 4)

-k 'number' : The support size in the superresolved grid. For very sharp
              images 4s+1 should be enough. (default 4s + 1)

-o 'filename' : Save the estimated PSF as an 8-bit image file. Values are
                rescaled to max_value be 255. Just for visualization purposes.

-p 'filename' : Input Noise Pattern PGM Image.  (default pattern_noise.pgm)

-d 'filename' : Save the input image into a PGM image file 'filename' with
                marks
                where the PSF pattern has been detected.

-t <0,1,2> : Choose the numerical algorithm for solving Least Squares.
                     0 - Least Squares
                     1 - Least Squares + thresholding  (default)
                     2 - Non-negative least Squares (slowest)


<input file> : Input 8/16 bits PGM Image

<output file> : Output PSF written into a TXT file as a k x k matrix of floats


Example
./psf_estim img_example.pgm psf.txt



Documentation
-------------
The following is an implementation of the Point Spread Function Estimation
Algorithm presented in:

[1] M. Delbracio, P. Muse, A. Almansa and J.M. Morel.
 The non-parametric sub-pixel local point spread function estimation is a well
 posed problem.  Submitted to the International Journal of Computer
 Vision (IJCV), November 2010.

and in more detail described on the portal IPOL www.ipol.im where there
 is much more information, including this code and an online demo version:

http://www.ipol.im/pub/algo/admm_non_blind_psf_estimation/


Please report bugs in psf_estim to <mdelbra@gmail.com>.

BUGS / HISTORY
--------------

Changes in Version 1.3
--> indentation fixed. TABs replaced by four spaces

Changes in Version 1.2
--> indentation fixed to four spaces
--> style code was changed according to reviewers from the IPOL journal

Changes in Version 1.1
--> bicubic parameter 'a' was wrongly set to 0.5 instead of -0.5.
--> thin-plate regularization parameter is from [0,Inf) instead of [0,1]
--> thin-plate regularization parameter is set to 10 (instead of 0.1)
--> order of parameters nx and ny in 'fftwf_plan_r2r_2d' where in reverse
    order. First should be ny then nx (as we work in Raw major)


Copyright and License
---------------------

"Point Spread Function Estimation from a Random Target"

Copyright 2010-2011 mauricio delbracio (mdelbra@gmail.com)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
