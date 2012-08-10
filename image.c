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
 * @file image.c
 * @brief library code for basic image processing.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <float.h>
#include "image.h"


/*Version 1.2 24 November 2011*/

/** @brief Error/Exit print a message and exit.
 *  @param msg
 */
static void error(char *msg)
{
    fprintf(stderr, "PSF_ESTIM Error: %s\n", msg);
    exit(EXIT_FAILURE);
}


/** @brief Euclidean distance between two points.
 *  @param x1
 *  @param y1
 *  @param x2
 *  @param y2
 *  @return distance between point (x1,y1) and (x2,y2)
 */
float dist_l2(float x1, float y1, float x2, float y2)
{
    return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}


/**
 * @brief Free memory used in ImageFloat 'i'
 * @param i
 */
void free_imageFloat(ImageFloat i)
{
    if (i == NULL || i->val == NULL)
    error("free_image_float: invalid input image.");
    free((void *) i->val);
    free((void *) i);
}

/**
 * @brief Create new ImageFloat of size 'nrow' x 'ncol'
 * @param ncol - number of columns
 * @param nrow - number of rows
 * @return created ImageFloat
 */
ImageFloat new_imageFloat(int ncol, int nrow)
{
    ImageFloat image;

    if (ncol == 0 || nrow == 0)
    error("new_image_float: invalid image size.");

    image = (ImageFloat) malloc(sizeof(struct imageFloatStruct));
    if (image == NULL)
    error("not enough memory.");
    image->val = (float *) calloc(ncol * nrow, sizeof(float));
    if (image->val == NULL)
    error("not enough memory.");

    image->ncol = ncol;
    image->nrow = nrow;

    return image;
}



/**
 * @brief Separable full convolution between an ImageFloat and a kernel
 * @param in - number of columns
 * @param xker - horizontal kernel (array of floats)
 * @param xsize - length of horizontal kernel
 * @param yker - vertical kernel (array of floats)
 * @param ysize - length of vertical kernel
 * @return convolved ImageFloat
 */
ImageFloat convol_sep2(ImageFloat in, float *xker, int xsize, float *yker,
                       int ysize)
{
    ImageFloat aux, out;
    int N, M, n, m, q;
    int h;
    float sum;

    /* get memory for images */
    N = in->ncol + xsize - 1;
    M = in->nrow + ysize - 1;
    aux = new_imageFloat(N, in->nrow);
    out = new_imageFloat(N, M);

    /*convolution boundry value: 0 */

    /* First sampling: x axis */
    for (n = 0; n < aux->ncol; n++)
    {
        for (m = 0; m < aux->nrow; m++)
        {
            sum = 0.0;
            for (q = 0; q < xsize; q++)
            {
                h = n - q;
                /* null boundary condition */
                if (h >= 0 && h < in->ncol)
                    sum += in->val[h + m * in->ncol] * xker[q];
            }
            aux->val[n + m * aux->ncol] = sum;
        }
    }

    /* Second sampling: y axis */
    for (m = 0; m < out->nrow; m++)
    {
        for (n = 0; n < out->ncol; n++)
        {
            sum = 0.0;
            for (q = 0; q < ysize; q++)
            {
                h = m - q;

                /* null boundary condition */
                if (h >= 0 && h < in->nrow)
                    sum += aux->val[n + h * aux->ncol] * yker[q];
            }
            out->val[n + m * out->ncol] = sum;
        }
    }

    /* free memory */ ;
    free_imageFloat(aux);

    return out;
}


/**
 * @brief Extract subimage from  integer pixels 'xmin' to 'xmax'
 *        and 'ymin' to 'ymax'
 * @param in - input ImageFloat
 * @param xmin
 * @param xmax
 * @param ymin
 * @param ymax
 * @return ImageFloat created from [xmin,xmax]X[ymin,ymax] region
 */
ImageFloat extract_window(ImageFloat in, int xmin,
                          int xmax, int ymin, int ymax)
{
    ImageFloat out;
    int i, j;

    out = new_imageFloat(xmax - xmin + 1, ymax - ymin + 1);

    for (i = 0; i < out->nrow; i++)
    for (j = 0; j < out->ncol; j++)
    {
        out->val[j + i * out->ncol] =
        in->val[j + xmin + (i + ymin) * in->ncol];
    }

    return out;
}

/**
 * @brief Extract subimage centered in a non-integer pixel (cx,cy)
 *        and square window of size 2*wsize + 1
 * @param in - input ImageFloat
 * @param wsize - window size (half side of the square)
 * @param cx - Center float x-coordinate
 * @param cy - Center float y-coordinate
 * @return ImageFloat square extracted from 'in' at (cx,cy) and side 2*wsize + 1
 */
ImageFloat extract_subpx_window(ImageFloat in, int wsize, float cx,
                                float cy)
{
    int cx_inf = (int) floor(cx);
    int cy_inf = (int) floor(cy);
    float dx = cx - cx_inf;
    float dy = cy - cy_inf;
    int xmin, xmax, ymin, ymax;
    ImageFloat out, aux;

    ImageFloat win;

    float vIx[2];
    float vIy[2];

    vIx[0] = dx;
    vIx[1] = 1 - dx;
    vIy[0] = dy;
    vIy[1] = 1 - dy;

    xmin = cx_inf - wsize - 1;
    xmax = cx_inf + wsize + 1;
    ymin = cy_inf - wsize - 1;
    ymax = cy_inf + wsize + 1;

    if (xmin < 0)
        xmin = 0;
    if (xmax >= in->ncol)
        xmax = in->ncol - 1;
    if (ymin < 0)
        ymin = 0;
    if (ymax >= in->nrow)
        ymax = in->nrow - 1;

    /*Extract the window of interest and a little more */
    win = extract_window(in, xmin, xmax, ymin, ymax);
    aux = convol_sep2(win, vIx, 2, vIy, 2);

    /*Extract the window of interest */
    out = extract_window(aux, 2, 2 * wsize + 2, 2, 2 * wsize + 2);

    free_imageFloat(win);
    free_imageFloat(aux);

    return out;
}


/**
 * @brief Horizontal image gradient calculated by finite differences
 * @param in - input ImageFloat
 * @return ImageFloat with the horizontal gradient
 */
ImageFloat gradx(ImageFloat in)
{
    int i, j;
    ImageFloat out;
    out = new_imageFloat(in->ncol, in->nrow);

    /* out-of-boundary calculus */
    for (i = 0; i < out->nrow; i++)
        for (j = 1; j < out->ncol - 1; j++)
            out->val[j + i * out->ncol] = - 0.5 *
                    in->val[j - 1 + i * out->ncol]
                    + 0.5 * in->val[j + 1 + i * out-> ncol];


    /* in-boundary calculus */
    j = 0;
    for (i = 0; i < out->nrow; i++)
        out->val[j + i * out->ncol] = - 0.5 * in->val[j + i * out->ncol]
            + 0.5 * in->val[j + 1 +i *out->ncol];

    j = out->ncol - 1;
    for (i = 0; i < out->nrow; i++)
        out->val[j + i * out->ncol] = - 0.5 * in->val[j - 1 + i * out->ncol]
            + 0.5 * in->val[j + i * out->ncol];

    return out;
}

/**
 * @brief Vertical image gradient calculated by finite differences
 * @param in - input ImageFloat
 * @return ImageFloat with the vertical gradient
 */
ImageFloat grady(ImageFloat in)
{
    int i, j;
    ImageFloat out;
    out = new_imageFloat(in->ncol, in->nrow);

    /* out-of-boundary calculus */
    for (i = 1; i < out->nrow - 1; i++)
        for (j = 0; j < out->ncol; j++)
            out->val[j + i * out->ncol] = - 0.5 *
                  in->val[j + (i - 1) * out->ncol]
                     + 0.5 * in->val[j + (i + 1)* out-> ncol];

    /* in-boundary calculus */
    i = 0;
    for (j = 0; j < out->ncol; j++)
        out->val[j + i * out->ncol] = - 0.5 * in->val[j + i * out->ncol]
            + 0.5 * in->val[j + (i + 1) * out->ncol];

    i = out->nrow - 1;
    for (j = 0; j < out->ncol; j++)
        out->val[j + i * out->ncol] = - 0.5 * in->val[j + (i - 1) * out->ncol]
            + 0.5 * in->val[j + i * out-> ncol];

    return out;
}


/**
 * @brief Bilinear image interpolation at coordinates given by 'X' and 'Y'.
 * @details It is assumed that the input image is sampled in a rectangular grid
 *          from (0,nc-1)X(0,nr-1). X and Y must be the same size and the
 *            output will be also the same size.
 * @param X - input ImageFloat with the horizontal coordinates where the
 *            interpolation is demanded
 * @param Y - input ImageFloat with the vertical coordinates where the
 *            interpolation is demanded
 * @param in - input ImageFloat
 * @return ImageFloat with the interpolated values at (X(i,j),Y(i,j))
 *         positions
 */
ImageFloat bilinear(ImageFloat X, ImageFloat Y, ImageFloat in)
{
    /* X and Y are the coordinates where the output image is going to be
     * interpolated. Both images should be of the same time and equal the
     * size of the output image. I assume that image in is sampled in a
     * rectangular regular grid from x=0:nc-1 and y=0:nr-1.
     */

    ImageFloat out;

    int nco = X->ncol;
    int nro = X->nrow;
    int nci = in->ncol;
    int nri = in->nrow;

    int xinf, yinf;
    float u, v, xp, yp, res;

    int i, j;

    /*Create Output image */
    out = new_imageFloat(nco, nro);


    for (j = 0; j < nco; j++)
    {
        for (i = 0; i < nro; i++)
        {
            xp = X->val[i * nco + j];
            yp = Y->val[i * nco + j];

            xinf = floor(xp);
            yinf = floor(yp);

            u = xp - (float) xinf;
            v = yp - (float) yinf;

            if (xinf > nci - 2 || xinf < 0 || yinf > nri - 2 || yinf < 0)
                res = 0;
            else
                res = (1 - u) * (1 - v) * in->val[yinf * nci + xinf] +
                    (1 - u) * v * in->val[(yinf + 1) * nci + xinf] +
                    u * (1 - v) * in->val[yinf * nci + xinf + 1] +
                    u * v * in->val[(yinf + 1) * nci + xinf + 1];

            out->val[i * nco + j] = res;
        }
    }

    return out;

}

/**
 * @brief Bicubic image interpolation at coordinates given by 'X' and 'Y'.
 * @details It is assumed that the input image is sampled in a rectangular grid
 *          from (0,nc-1)X(0,nr-1). X and Y must be the same size and
 *          the output will be also the same size. Where there is not enough
 *          information to interpolate (i.e. in the boundary) the function
 *          returns 0
 *
 * @param X - input ImageFloat with the horizontal coordinates where the
 *            interpolation is demanded
 * @param Y - input ImageFloat with the vertical coordinates where the
 *        interpolation is demanded
 * @param in - input ImageFloat
 * @param a -  Bicubic inerpolator parameter (tipically -0.5)
 * @return ImageFloat with the interpolated values at (X(i,j),Y(i,j))
 *         positions
 */
ImageFloat bicubic(ImageFloat X, ImageFloat Y, ImageFloat in, float a)
{
    /* X and Y are the coordinates where the output image is going to be
     * interpolated. Both images should be of the same time and equal the
     * size of the output image. I assume that image in is sampled in a
     * rectangular regular grid from x=0:nc-1 and y=0:nr-1. Where there is
     * not enough information to interpolate (i.e. in the boundary) the
     * function returns 0.
     */

    ImageFloat out;

    int nco = X->ncol;
    int nro = X->nrow;
    int nci = in->ncol;
    int nri = in->nrow;

    float cx[4], cy[4];

    int xinf, yinf;
    float xp, yp, res, t, s, at, as, t2, s2;

    int i, j;

    /*Create Output image */
    out = new_imageFloat(nco, nro);

    /*----Bicubic interpolation - Central Loop----*/

    for (j = 0; j < nco; j++)
    {
        for (i = 0; i < nro; i++)
        {
            xp = X->val[i * nco + j];
            yp = Y->val[i * nco + j];

            xinf = floor(xp);
            yinf = floor(yp);

            t = xp - (float) xinf;
            t2 = t * t;
            at = a * t;
            cx[0] = a * t2 * (1.0 - t);
            cx[1] = (2.0 * a + 3.0 - (a + 2.0) * t) * t2 - at;
            cx[2] = ((a + 2.0) * t - a - 3.0) * t2 + 1.0;
            cx[3] = a * (t - 2.0) * t2 + at;
            s = yp - (float) yinf;
            s2 = s * s;
            as = a * s;
            cy[0] = a * s2 * (1.0 - s);
            cy[1] = (2.0 * a + 3.0 - (a + 2.0) * s) * s2 - as;
            cy[2] = ((a + 2.0) * s - a - 3.0) * s2 + 1.0;
            cy[3] = a * (s - 2.0) * s2 + as;

            /*Put 0 in the border.... */
            if (xinf > nci - 3 || xinf < 1 || yinf > nri - 3 || yinf < 1)
                res = 0;
            else
            /* Re-check if this is not separable...or something */
                res = cy[0] * (cx[0] * in->val[(yinf + 2) * nci + xinf + 2] +
                               cx[1] * in->val[(yinf + 2) * nci + xinf + 1] +
                               cx[2] * in->val[(yinf + 2) * nci + xinf] +
                               cx[3] * in->val[(yinf + 2) * nci + xinf - 1]
                               ) +
                    cy[1] * (cx[0] * in->val[(yinf + 1) * nci + xinf + 2] +
                             cx[1] * in->val[(yinf + 1) * nci + xinf + 1] +
                             cx[2] * in->val[(yinf + 1) * nci + xinf] +
                             cx[3] * in->val[(yinf + 1) * nci + xinf - 1]
                             ) +
                    cy[2] * (cx[0] * in->val[yinf * nci + xinf + 2] +
                             cx[1] * in->val[yinf * nci + xinf + 1] +
                             cx[2] * in->val[yinf * nci + xinf] +
                             cx[3] * in->val[yinf * nci + xinf - 1]
                             ) +
                    cy[3] * (cx[0] * in->val[(yinf - 1) * nci + xinf + 2] +
                             cx[1] * in->val[(yinf - 1) * nci + xinf + 1] +
                             cx[2] * in->val[(yinf - 1) * nci + xinf] +
                             cx[3] * in->val[(yinf - 1) * nci + xinf - 1]
                             );
            out->val[i * nco + j] = res;
        }
    }

    return out;
}

/**
 * @brief Calculate the Power (mean of the square image) of a subimage
 *        (xmin,xmax) X (ymin,ymax)
 * @param in - input ImageFloat
 * @param xmin
 * @param xmax
 * @param ymin
 * @param ymax
 * @return mean of the square selected subimage
 */
float power_window(ImageFloat in, int xmin, int xmax, int ymin, int ymax)
{
    float out = 0;
    int i, j;

    for (i = ymin; i < ymax + 1; i++)
        for (j = xmin; j < xmax + 1; j++)
        {
            out += (in->val[j + i * in->ncol]) * (in->val[j + i * in->ncol]);
        }

    out = out / ((xmax - xmin + 1) * (ymax - ymin + 1));

    return out;
}

/**
 * @brief Calculate the mean value of a subimage (xmin,xmax) X (ymin,ymax)
 * @param in - input ImageFloat
 * @param xmin
 * @param xmax
 * @param ymin
 * @param ymax
 * @return mean value of the selected subimage
 */
float mean_window(ImageFloat in, int xmin, int xmax, int ymin, int ymax)
{
    float out = 0;
    int i, j;

    for (i = ymin; i < ymax + 1; i++)
        for (j = xmin; j < xmax + 1; j++)
        {
            out += (in->val[j + i * in->ncol]);
        }

    out = out / ((xmax - xmin + 1) * (ymax - ymin + 1));

    return out;
}


/**
 * @brief Calculate the mean value of a subimage at center in non-integer pixel
 *        (cx,cy) and square window of size 2*wsize + 1
 * @param in - input ImageFloat
 * @param wsize - window size (half side of the square)
 * @param cx - Center float x-coordinate
 * @param cy - Center float y-coordinate
 * @return mean value of the selected subimage
 */
float mean_subpx_window(ImageFloat in, int wsize, float cx, float cy)
{
    int cx_inf = (int) floor(cx);
    int cy_inf = (int) floor(cy);
    float dx = cx - cx_inf;
    float dy = cy - cy_inf;
    int xmin, xmax, ymin, ymax;
    ImageFloat aux;
    float out;

    float vIx[2], vIy[2];

    vIx[0] = dx;
    vIx[1] = 1 - dx;
    vIy[0] = dy;
    vIy[1] = 1 - dy;

    aux = convol_sep2(in, vIx, 2, vIy, 2);

    xmin = cx_inf - wsize + 1;
    xmax = cx_inf + wsize + 1;
    ymin = cy_inf - wsize + 1;
    ymax = cy_inf + wsize + 1;

    if (xmin < 0)
        xmin = 0;
    if (xmax >= aux->ncol)
        xmax = aux->ncol - 1;
    if (ymin < 0)
        ymin = 0;
    if (ymax >= aux->nrow)
        ymax = aux->nrow - 1;

    /*printf("%d %d %d %d\n",xmin,xmax,ymin,ymax); */
    out = mean_window(aux, xmin, xmax, ymin, ymax);

    free_imageFloat(aux);

    return out;
}

/**
 * @brief Calculate the Power (mean value of the square) of a subimage at
 *        center in non-integer pixel (cx,cy) and square window of size
 *        2*wsize + 1
 * @param in - input ImageFloat
 * @param wsize - window size (half side of the square)
 * @param cx - Center float x-coordinate
 * @param cy - Center float y-coordinate
 * @return mean value of the square of selected subimage (power)
 */
float power_subpx_window(ImageFloat in, int wsize, float cx, float cy)
{

    int cx_inf = (int) floor(cx);
    int cy_inf = (int) floor(cy);
    float dx = cx - cx_inf;
    float dy = cy - cy_inf;
    int xmin, xmax, ymin, ymax;
    ImageFloat aux;
    float out;

    float vIx[2], vIy[2];

    vIx[0] = dx;
    vIx[1] = 1 - dx;
    vIy[0] = dy;
    vIy[1] = 1 - dy;

    aux = convol_sep2(in, vIx, 2, vIy, 2);

    xmin = cx_inf - wsize + 1;
    xmax = cx_inf + wsize + 1;
    ymin = cy_inf - wsize + 1;
    ymax = cy_inf + wsize + 1;

    if (xmin < 0)
        xmin = 0;
    if (xmax >= aux->ncol)
        xmax = aux->ncol - 1;
    if (ymin < 0)
        ymin = 0;
    if (ymax >= aux->nrow)
        ymax = aux->nrow - 1;

    out = power_window(aux, xmin, xmax, ymin, ymax);

    free_imageFloat(aux);

    return out;
}

/**
 * @brief Compute the DCT Transform of ImageFloat 'in'
 * @param in - input ImageFloat
 * @return ImageFloat computed DCT image
 */
ImageFloat compute_dct_image(ImageFloat in)
{

    ImageFloat out;
    fftwf_plan fw_plan;
    int nx, ny;
    nx = in->ncol;
    ny = in->nrow;

    out = new_imageFloat(nx, ny);

    /*Be careful! the order of the parameters: ny then nx
     "The multi-dimensional arrays passed to fftw_plan_dft etcetera are
     expected to be stored as a single contiguous block in row-major order
     (sometimes called “C order”).
     Basically, this means that as you step through adjacent
     memory locations, the first dimension's index varies most slowly and the
     last dimension's index varies most quickly."*/
    fw_plan = fftwf_plan_r2r_2d(ny, nx, in->val, out->val, FFTW_REDFT10,
                                FFTW_REDFT10, FFTW_ESTIMATE);
    fftwf_execute(fw_plan);

    /* Do the cleaning */
    fftwf_destroy_plan(fw_plan);

    return out;
}

/**
 * @brief Compute  the DCT Inverse Transform of ImageFloat 'in'
 * @param in - input ImageFloat
 * @return ImageFloat computed DCT image
 */
ImageFloat compute_idct_image(ImageFloat in)
{

    ImageFloat out;
    fftwf_plan bw_plan;
    int nx, ny;
    nx = in->ncol;
    ny = in->nrow;

    out = new_imageFloat(nx, ny);

    /*Be careful! the order of the parameters: ny then nx
     "The multi-dimensional arrays passed to fftw_plan_dft etcetera are
     expected to be stored as a single contiguous block in row-major order
     (sometimes called “C order”).
     Basically, this means that as you step through adjacent
     memory locations, the first dimension's index varies most slowly and the
     last dimension's index varies most quickly."*/

    bw_plan = fftwf_plan_r2r_2d(ny, nx, in->val, out->val, FFTW_REDFT01,
                                FFTW_REDFT01, FFTW_ESTIMATE);

    fftwf_execute(bw_plan);

    /* Do the cleaning */
    fftwf_destroy_plan(bw_plan);

    return out;
}
