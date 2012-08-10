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
 * @file image.h
 * @brief library header for basic image processing.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */

#ifndef IMAGE_H_
#define IMAGE_H_




/*---------------------------------------------------------------------------*/
/** @brief double image val type

    The pixel value at (x,y) is accessed by:

    image->val[ x + y * image->ncol ]

    with x and y integer.
 */
typedef struct imageFloatStruct
{
    float *val;
    int ncol, nrow;
} *ImageFloat;

void free_imageFloat(ImageFloat i);
ImageFloat new_imageFloat(int ncol, int nrow);


/*---------------------------------------------------------------------------*/




float dist_l2(float x1, float y1, float x2, float y2);

ImageFloat convol_sep2(ImageFloat in, float *xker, int xsize, float *yker,
                       int ysize);


ImageFloat extract_window(ImageFloat in, int xmin, int xmax,
                          int ymin, int ymax);

ImageFloat extract_subpx_window(ImageFloat in, int wsize, float cx,
                                float cy);

ImageFloat gradx(ImageFloat in);

ImageFloat grady(ImageFloat in);

ImageFloat bilinear(ImageFloat X, ImageFloat Y, ImageFloat in);

ImageFloat bicubic(ImageFloat X, ImageFloat Y, ImageFloat in, float a);

float mean_window(ImageFloat in, int xmin, int xmax, int ymin, int ymax);

float power_window(ImageFloat in, int xmin, int xmax, int ymin, int ymax);


float mean_subpx_window(ImageFloat in, int wsize, float cx, float cy);
float power_subpx_window(ImageFloat in, int wsize, float cx, float cy);

ImageFloat compute_dct_image(ImageFloat img);
ImageFloat compute_idct_image(ImageFloat img);


#endif          /* IMAGE_H_ */
