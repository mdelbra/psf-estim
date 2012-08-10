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
 * @file thin_plates.h
 * @brief library header to estimate/evaluate thin plates splines.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#ifndef THIN_PLATES_H_
#define THIN_PLATES_H_

#define MAX_NC 1000    /**< Maximum number of allowed points */
#define D  3/**< Dimensions: (x,y,1) in homogeneaus coordinates */
#define EPS_ZERO 0.00000000001 /**< Very small value for zero tolerance */
#define BUFFER_SIZE 113337 /**< Buffer size defined empirically to consider
                 * up to 1000 points */

/** Thin Plate Structure */
typedef struct thinPlateStruct {
    int nc;             /**< number of centers */
    float *xc;          /** x-coordinate of centers */
    float *yc;          /** y-coordinate of centers */
    float *coef_x;      /** non-affine coeficient for x-output */
    float *coef_y;      /** non-affine coeficient for y-output */
    float affine[6];    /** affine coeficients
                            0,1,2 for x_output part
                            3,4,5 for y_output part */
    float lambda;   /** regularization parameter - just for the record */
} *ThinPlate;


ThinPlate calculate_thinPlate(float *Pin, float *Pout, int k,
                              float lambda);

int evaluate_thinPlate(ThinPlate tp, float *Pin, float *Pout, int np);

void free_thinPlate(ThinPlate tp);

#endif                /* THIN_PLATES_H_ */
