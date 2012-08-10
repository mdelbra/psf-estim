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
 * @file detect_pattern.h
 * @brief  Module header for detecting the pattern in a digital image
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#ifndef DETECT_PATTERN_H_
#define DETECT_PATTERN_H_


#define EXIT_SEVERAL_PATTERNS_DETECTED -100
#define EXIT_NO_PATTERN_DETECTED -101


/* The analytic pattern is rasterized 8x8 per black/white value */
#define UP_RES 8.0

/* The number of black and white pattern values per block (checkerboard)*/
#define PATTERN_BLOCK_SIZE 32.0


float *detect_pattern(ImageFloat in);

ImageFloat draw_detected_corners_image_maxval(ImageFloat in,
                                              float *checkerboard);
ImageFloat draw_detected_corners_image_minval(ImageFloat in,
                                              float *checkerboard);

float *pattern_Xpoints(void);
void pattern_whiteSquares(float *p);
void pattern_blackSquares(float *p);
void pattern_center(float *p);
void pattern_top_center(float *p);

int roundfi(float x);

ImageFloat pattern_to_pattern_image(float* pattern, int pat_nx, int pat_ny);

#endif          /* DETECT_PATTERN_H_ */
