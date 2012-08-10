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


/**
 * @file nnls.h
 * @brief library header numerical algorithms for solving least squares
 *        and non-negative least squares
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 * @date Oct 11, 2010
 */


#ifndef NNLS_H_
#define NNLS_H_


float *solve_nnls(float *A, float *b, int m, int n);
float *solve_ls_th(float *A, float *b, int m, int n, float th);
float *solve_ls(float *A, float *b, int m, int n);

float *solve_nnlsd(float *A, float *b, int m, int n);
float *solve_lsd_th(float *A, float *b, int m, int n, float th);
float *solve_lsd(float *A, float *b, int m, int n);

#endif                /* NNLS_H_ */
