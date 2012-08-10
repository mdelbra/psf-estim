/*----------------------------------------------------------------------------

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

 ----------------------------------------------------------------------------*/


/*
 * psf_estim.h
 *
 *  Created on: Oct 22, 2010
 *      Author: mdelbra
 */

#ifndef PSF_ESTIM_H_
#define PSF_ESTIM_H_


void write_ascii_matrix(float *M, int ncol, int nrow, char *name);

float *psf_estim(float *img, int nx, int ny,
         float *pattern, int pat_nx, int pat_ny,
         int s, int psf_nrow, int psf_ncol, int solver, char* detected_pgm);

#endif                /* PSF_ESTIM_H_ */
