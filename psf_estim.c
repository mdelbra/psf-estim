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
 * @file psf_estim.c
 * @brief library code to psf estimation.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


/** @mainpage Non-parametric sub-pixel local point spread function estimation
 *
 * The following is an implementation of the Point Spread Function Estimation
 * Algorithm presented in
 *
 * \li M. Delbracio, P. Musé, A. Almansa and J.M. Morel.
 * The non-parametric sub-pixel local point spread function estimation is
 * a well posed problem. Submitted to the International Journal
 * of Computer Vision (IJCV), November 2010.
 *
 * and in more detail described on the portal IPOL www.ipol.im where
 * there is much more information, including this code and an
 * online demo version:
 *
 * \li http://www.ipol.im/pub/algo/admm_non_blind_psf_estimation/
 *
 * The source code consists of:
 *
 *    \li    detect_pattern.c
 *    \li    detect_pattern.h
 *    \li    image.c
 *    \li    image.h
 *    \li    io_pgm.c
 *    \li    io_pgm.h
 *    \li    lsd.c
 *    \li    lsd.h
 *    \li    nnls.c
 *    \li    nnls.h
 *    \li    psf_estim.c
 *    \li    psf_estim.h
 *    \li    psf_estim_main.c
 *    \li    psf_estim_main.h
 *    \li    thin_plates.c
 *    \li    thin_plates.h
 *
 *
 * HISTORY:
 * - Version 1.3 - 26 December 2011
 * - Version 1.2 - 24 Novemeber 2011
 * - Version 1.1 - 02 July 2011
 *
 *
 * @author mauricio delbracio (mdelbra@gmail.com)
 * @date jul 2011
 */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "image.h"
#include "detect_pattern.h"
#include "thin_plates.h"
#include "nnls.h"
#include "io_pgm.h"


/* For debugging purposes you can save all intermediate images. By default
 this is off.
 #DEFINE SAVE_INTERMEDIATE_IMGS
 */

/** @brief Error/Exit print a message and exit.
 *  @param msg
 */
void error (char *msg)
{
    fprintf (stderr, "psf_estim Error: %s\n", msg);
    exit (EXIT_FAILURE);
}


/**
 *  @brief Write an ImageFloat to an ASCII file
 *  @param image - image to write
 *  @param name  - filename
 */
void write_ascii_imageFloat (ImageFloat image, char *name)
{
    FILE *f;
    int x, y, n;

    /* open file */
    f = fopen (name, "w");
    if (f == NULL)
    {
        error ("Can't open output file.");
    }

    /* write data */
    for (y = 0; y < image->nrow; y++)
    {
        for (x = 0; x < image->ncol; x++, n++)
            fprintf (f, "%f ", image->val[x + y * image->ncol]);
        fprintf (f, "\n");
    }

    /* close file */
    fclose (f);
}


/**
 *  @brief Write an Array of floats (Matrix)  to an ASCII file
 *  @param M - array of floats to write
 *  @param ncol - number of columns of M
 *  @param nrow - number of rows of M
 *  @param name  - filename
 */
void write_ascii_matrix(float *M, int ncol, int nrow, char *name)
{
    FILE *f;
    int x, y, n;

    /* open file */
    f = fopen (name, "w");
    if (f == NULL)
        error ("Can't open output file.");

    /* write header */

    /* write data */
    for (y = 0; y < nrow; y++)
    {
        for (x = 0; x < ncol; x++, n++)
            fprintf (f, "%f ", M[y + x * nrow]);

        fprintf (f, "\n");
    }

    /* close file */
    fclose (f);
}


/**
 *  @brief Image Normalization
 */
ImageFloat image_normalization (ImageFloat in, ThinPlate tp)
{
    /*Image Normalization - BLACK */
    /*----------------------------*/

    /* Coordinates of the center white squares in Sharp and Blur images */
    /* There are 12 black and 12 white squares */
    float cwhiteS[24], cwhiteB[24];
    float cblackS[24], cblackB[24];

    /*To generate the matrix with 6 terms: x^2 y^2 xy x y 1 */
    float CblackM[12 * 6], CwhiteM[12 * 6];
    float xc, yc;

    /* Gray value (mean)of the square black and white regions */
    float gblack[12], gwhite[12], *pol_white, *pol_black, xW, yW, xB, yB, x, y;


    int np = 12, wsize = 2, i, j;

    int nc, nr;

    ImageFloat img_black, img_white, out;

    pattern_whiteSquares (cwhiteS);
    pattern_blackSquares (cblackS);

    evaluate_thinPlate (tp, cwhiteS, cwhiteB, np);
    evaluate_thinPlate (tp, cblackS, cblackB, np);


    /* Arbitrarily we set the (0,0) of the coordinate system in the center
     * of the cwhiteB points. If I do not do this and I put the (0,0) in
     * one particular point the System will adjust better to some regions
     * Also, the coordinates values are rounded in order to access the pixel
     * value of the same
     */

    for (i = 0, xc = 0, yc = 0; i < np; i++)
    {
        xc += cwhiteB[2 * i];
        yc += cwhiteB[2 * i + 1];
    }

    xc /= np;
    yc /= np;


    for (i = 0; i < np; i++)
    {
        xW = cwhiteB[2 * i] - xc;
        yW = cwhiteB[2 * i + 1] - yc;
        xB = cblackB[2 * i] - xc;
        yB = cblackB[2 * i + 1] - yc;

        /*gblackM and gwhiteM are matrix with coordinate locations */
        CblackM[i] = xB * yB;
        CblackM[i + np] = xB;
        CblackM[i + 2 * np] = yB;
        CblackM[i + 3 * np] = 1;
        CblackM[i + 4 * np] = xB * xB;
        CblackM[i + 5 * np] = yB * yB;

        CwhiteM[i] = xW * yW;
        CwhiteM[i + np] = xW;
        CwhiteM[i + 2 * np] = yW;
        CwhiteM[i + 3 * np] = 1;
        CwhiteM[i + 4 * np] = xW * xW;
        CwhiteM[i + 5 * np] = yW * yW;


        /*Compute the mean gray values in a window center in the square */
        gblack[i] =
        mean_subpx_window (in, wsize, cblackB[2 * i], cblackB[2 * i + 1]);
        gwhite[i] =
        mean_subpx_window (in, wsize, cwhiteB[2 * i], cwhiteB[2 * i + 1]);

    }

    /*Solving the Linear System CblackM * x = gB and CwhiteM * x = gW */
    /* pol_black =  inv(CblackM)*gblack   --  12 x 6 */
    pol_black = solve_ls(CblackM,gblack,12,6);

    /* gwhite =  inv(CwhiteM)*gwhite   --  12 x 6 */
    pol_white = solve_ls(CwhiteM,gwhite,12,6);


    /* Generating Black and White illumination images */

    nc = in->ncol;
    nr = in->nrow;

    /*Strictly this is not necessary as I can do directly the normalization
     * without generating the black and white images
     */

    img_black = new_imageFloat (nc, nr);
    img_white = new_imageFloat (nc, nr);
    out = new_imageFloat (nc, nr);


    for (i = 0; i < nr; i++)
        for (j = 0; j < nc; j++)
        {
            x = j - xc;
            y = i - yc;
            img_white->val[j + i * nc] = x * y * pol_white[0]
                                        + x * pol_white[1]
                                        + y * pol_white[2]
                                        + pol_white[3]
                                        + x * x * pol_white[4]
                                        + y * y * pol_white[5];
            img_black->val[j + i * nc] = x * y * pol_black[0]
                                        + x * pol_black[1]
                                        + y * pol_black[2]
                                        + pol_black[3]
                                        + x * x * pol_black[4]
                                        + y * y * pol_black[5];

            out->val[j + i * nc] =  (in->val[j + i * nc]
                                     - img_black->val[j + i * nc])
                                        / (img_white->val[j + i * nc]
                                     - img_black->val[j + i * nc]);
        }



    free_imageFloat (img_white);
    free_imageFloat (img_black);
    free((void *) pol_white);
    free((void *)    pol_black);

    return out;

}


/**
 *  @brief Camera Response Function Correction
 */
ImageFloat crf_correction (ImageFloat in, ThinPlate tp)
{
    /*Iterative adjustment to generate an image of mean(noise_part) = 0.5
     * This module supposes that the image is normalized between 0 and 1, where
     * 0 is the black value and 1 the white value.
     * */

    float mean, power, a;
    int i;
    float center_sharp[2], center_blur[2];
    float top_sharp[2], top_blur[2];
    int wsize;
    ImageFloat out;

    out = new_imageFloat (in->ncol, in->nrow);

    /*Calculating the center of the pattern in the blur image */
    pattern_center (center_sharp);
    evaluate_thinPlate (tp, center_sharp, center_blur, 1);

    pattern_top_center(top_sharp);
    evaluate_thinPlate (tp, top_sharp, top_blur, 1);

    /* approx 0.85 half the size of the noise region in the blur image*/
    wsize = (int)(0.85 * dist_l2(top_blur[0],top_blur[1],
                                 center_blur[0],center_blur[1]));


    /*A parabolic function is estimated as:
     * y = ax^2 + (1-a)x
     *  This parabola maps 0 -> 0, 1-> 1
     *  a is set in order that yM -> 0.5
     *  so a = (0.5-mean(x))/(mean(x^2) - mean(x))
     */

    mean = mean_subpx_window (in, wsize, center_blur[0], center_blur[1]);
    power = power_subpx_window (in, wsize, center_blur[0], center_blur[1]);
    a = (0.5 - mean) / (power - mean);

    for (i = 0; i < in->ncol * in->nrow; i++)
        out->val[i] = a * in->val[i] * in->val[i] + (1 - a) * in->val[i];

    return out;
}


/**
 *  @brief Low-Pass filter Image (with DCT Transform)
 */
ImageFloat lpf_image_dct (ImageFloat in, int fcx, int fcy)
{
    ImageFloat data_dct, out;
    int i, j;
    int nx = in->ncol;
    int ny = in->nrow;
    float k;

    data_dct = compute_dct_image (in);


    for (i = fcy; i < ny; i++)
        for (j = 0; j < nx; j++)
        {
            data_dct->val[i * nx + j] = 0;
        }


    for (i = 0; i < ny; i++)
        for (j = fcx; j < nx; j++)
        {
            data_dct->val[i * nx + j] = 0;
        }


    out = compute_idct_image (data_dct);

    /*Normalizing image because DCT introduces
     * a constant factor 4*nx*ny
     */
    k = 4 * nx * ny;
    for (i = 0; i < nx * ny; i++)
        out->val[i] = out->val[i] / k;

    /*
     * cleanup
     */
    free_imageFloat (data_dct);

    return out;

}

/**
 *  @brief Generate the Linear System Ax = b where x is the PSF to find.
 */
int make_Ab (ImageFloat imgC, ImageFloat imgW, ImageFloat imgMask,
             int q, int p, int s, float **A, float *b[], int *ncol, int *nrow)
{
    /*Generate a SsU from U and s. */
    int max_pq, u, v, i, j, mkMs, nkMs, mc, nc, r;

    float *kerode;
    ImageFloat mask;
    /* generate SsU matrix */
    /* p,q kernel size */

    /*size of the final matrix A: (mkMs*nKMs)x(p*q) */
    mkMs = (imgC->nrow + p - 2) / s + 1;
    nkMs = (imgC->ncol + q - 2) / s + 1;

    *A = (float *)calloc (mkMs * nkMs * p * q, sizeof (float));/*ini. to zero*/
    *b = (float *)calloc (mkMs * nkMs, sizeof (float));/*ini. to zero */


    /*Erode the Mask in order take into account the boundary
     * problems
     */
    /*square element of side r */
    max_pq = (p > q) ? p : q;
    r = (max_pq - 1) / (2 * s);
    kerode = (float *) malloc ((2 * r + 1) * sizeof (float));/*ini. to zero */

    for (i = 0; i < 2 * r + 1; i++)
        kerode[i] = 1;

    mask = convol_sep2 (imgMask, kerode, 2 * r + 1, kerode, 2 * r + 1);

    for (i = 0; i < imgMask->nrow; i++)
        for (j = 0; j < imgMask->ncol; j++)
            imgMask->val[i * imgMask->ncol + j] =
            (mask->val[(i + r) * mask->ncol + j + r] >= 4 * r * r + 1) ? 1 : 0;


    /*Filling in (*A) */
    for (i = (p - 1) / 2; i < (imgC->nrow) + (p - 1) / 2; i = i + s)
    {
        for (j = (q - 1) / 2; j < (imgC->ncol) + (q - 1) / 2; j = j + s)
        {
            if (imgMask->val[(j - (q - 1) / 2) / s +
                             imgMask->ncol * (i - (p - 1) / 2) / s])
            {
                for (u = 0; u < p; u++)
                {
                    for (v = 0; v < q; v++)
                    {
                        if ((i - u >= 0) && (j - v >= 0)
                            && (i - u < imgC->nrow)
                            && (j - v < imgC->ncol))
                            (*A)[mkMs * nkMs * (u + v * p)
                                 + (i / s + j / s * mkMs)]
                                = imgC->val[j - v + imgC->ncol * (i - u)];
                                  /*Save by cols*/
                    }
                }
            }
        }
    }

    /* Put the image in he center (just to have a centered kenel) */
    mc = (p - 1) / (2 * s) + 1;
    nc = (q - 1) / (2 * s) + 1;


    for (i = 0; i < imgW->nrow; i++)
        for (j = 0; j < imgW->ncol; j++)
        {
            /*Check if it is in the mask... */
            if (imgMask->val[i * imgMask->ncol + j])
                (*b)[(j + nc - 1) * mkMs + i + mc - 1] =
                    imgW->val[i * imgW->ncol + j];

        }


    *ncol = p * q;
    *nrow = mkMs * nkMs;

    free_imageFloat (imgMask);
    free_imageFloat (mask);
    free ((void *) kerode);

    return EXIT_SUCCESS;


}


/**
 *  @brief Extract noise Region from the observed image. Only extracts
 the random noise region for psf estimation (without the suroundings).
 */
ImageFloat extract_noise_region (ImageFloat in, ThinPlate tpI,
                                 float *xmin, float *xmax,
                                 float *ymin, float *ymax, ImageFloat * mask)
{

    ImageFloat imgW, imgMask;

    int i, j;
    float pwind_sharp[8], pwind_blur[8];

    float x, y, l1, l2, l3, xA, xB, xC, xD, yA, yB, yC, yD, detABC, detACD;
    char t1, t2;

    /*A*/
    pwind_sharp[0] = PATTERN_BLOCK_SIZE * UP_RES;
    pwind_sharp[1] = PATTERN_BLOCK_SIZE * UP_RES;

    /*B*/
    pwind_sharp[2] = PATTERN_BLOCK_SIZE * UP_RES;
    pwind_sharp[3] = PATTERN_BLOCK_SIZE * UP_RES * 9;

    /*C*/
    pwind_sharp[4] = PATTERN_BLOCK_SIZE * UP_RES * 9;
    pwind_sharp[5] = PATTERN_BLOCK_SIZE * UP_RES * 9;

    /*D*/
    pwind_sharp[6] = PATTERN_BLOCK_SIZE * UP_RES * 9;
    pwind_sharp[7] = PATTERN_BLOCK_SIZE * UP_RES;

    evaluate_thinPlate (tpI, pwind_sharp, pwind_blur, 4);

    /*Extract the smallest rectangle (parallel to x and y axis)
     * that includes the window noise region.
     */
    for (i = 1, *xmin = pwind_blur[0]; i < 4; i++)
        *xmin = (*xmin < pwind_blur[2 * i]) ? *xmin : pwind_blur[2 * i];

    for (i = 1, *xmax = pwind_blur[0]; i < 4; i++)
        *xmax = (*xmax > pwind_blur[2 * i]) ? *xmax : pwind_blur[2 * i];

    for (i = 1, *ymin = pwind_blur[1]; i < 4; i++)
        *ymin = (*ymin < pwind_blur[2 * i + 1]) ? *ymin : pwind_blur[2*i + 1];

    for (i = 1, *ymax = pwind_blur[1]; i < 4; i++)
        *ymax = (*ymax > pwind_blur[2 * i + 1]) ? *ymax : pwind_blur[2*i + 1];

    *xmin = floor (*xmin);
    *xmax = ceil (*xmax);
    *ymin = floor (*ymin);
    *ymax = ceil (*ymax);

    imgW =
        extract_window(in, (int) *xmin, (int) *xmax, (int) *ymin, (int) *ymax);


    imgMask = new_imageFloat (imgW->ncol, imgW->nrow);

    /* The region is the interior of a quadrilateral defined by
     * the four segments: AB, BC, CD, DA
     * http://en.wikipedia.org/wiki/Barycentric_coordinates_(mathematics)
     * */

    xA = pwind_blur[0];
    yA = pwind_blur[1];
    xB = pwind_blur[2];
    yB = pwind_blur[3];
    xC = pwind_blur[4];
    yC = pwind_blur[5];
    xD = pwind_blur[6];
    yD = pwind_blur[7];

    detABC = (xA - xC) * (yB - yC) - (xB - xC) * (yA - yC);
    detACD = (xA - xD) * (yC - yD) - (xC - xD) * (yA - yD);

    for (i = 0; i < imgW->nrow; i++)
        for (j = 0; j < imgW->ncol; j++)
        {
            x = *xmin + j;
            y = *ymin + i;
            l2 = ((yC - yA) * (x - xC) + (xA - xC) * (y - yC)) / detABC;
            l1 = ((yB - yC) * (x - xC) + (xC - xB) * (y - yC)) / detABC;
            l2 = ((yC - yA) * (x - xC) + (xA - xC) * (y - yC)) / detABC;
            l3 = 1 - l1 - l2;

            t1 = (0 <= l1) && (l1 <= 1) &&
            (0 <= l2) && (l2 <= 1) && (0 <= l3) && (l3 <= 1);

            /*If not in t1 check for t2 =  ACD */
            if (!t1)
            {
                l1 = ((yC - yD) * (x - xD) + (xD - xC) * (y - yD)) / detACD;
                l2 = ((yD - yA) * (x - xD) + (xA - xD) * (y - yD)) / detACD;
                l3 = 1 - l1 - l2;

                t2 = (0 <= l1) && (l1 <= 1) && (0 <= l2) && (l2 <= 1)
                        && (0 <= l3) && (l3 <= 1);
                if (t2)
                    imgMask->val[i * imgMask->ncol + j] = 1;
                else
                    imgMask->val[i * imgMask->ncol + j] = 0;
            }
            else
                imgMask->val[i * imgMask->ncol + j] = 1;
        }


    *mask = imgMask;

    return imgW;

}


/**
 *  @brief Extract pattern region from the observed image
 (noise part + surrunding marks, black/white squares)
 */
ImageFloat extract_pattern_region (ImageFloat in, float *p, int num_points)
{
    float xmin, xmax, ymin, ymax;
    int i;
    float e;
    ImageFloat imgT;

    /*Get the minimum rectangular window that covers all detected points */
    xmin = p[0];
    xmax = p[0];
    ymin = p[1];
    ymax = p[1];

    for (i = 1; i < num_points; i++)
    {
        xmin = (p[2 * i] < xmin) ? p[2 * i] : xmin;
        xmax = (p[2 * i] > xmax) ? p[2 * i] : xmax;
        ymin = (p[2 * i + 1] < ymin) ? p[2 * i + 1] : ymin;
        ymax = (p[2 * i + 1] > ymax) ? p[2 * i + 1] : ymax;
    }
    /*The rectangular window is enlarged a little to cover the whole pattern*/

    e = 0.1;

    xmin = floor (xmin - (xmax - xmin) * e);
    xmax = ceil (xmax + (xmax - xmin) * e);
    ymin = floor (ymin - (ymax - ymin) * e);
    ymax = ceil (ymax + (ymax - ymin) * e);

    /*Check values are in bound */
    xmin = xmin >= 0 ? xmin : 0;
    xmax = xmax < in->ncol ? xmax : in->ncol - 1;
    ymin = ymin >= 0 ? ymin : 0;
    ymax = ymax < in->nrow ? ymax : in->nrow - 1;

    imgT = extract_window (in, (int) xmin, (int) xmax, (int) ymin, (int) ymax);

    /*update checkpoints to the new references */

    for (i = 0; i < num_points; i++)
    {
        p[2 * i] -= xmin;
        p[2 * i + 1] -= ymin;

    }

    return imgT;

}



/**
 *  @brief PSF Estimation (Main Function)
 */
float * psf_estim (float *img, int nx, int ny,
                   float *pattern, int pat_nx, int pat_ny,
                   int s, int psf_nrow, int psf_ncol, int solver,
                   char *detected_ppm)
{
    ImageFloat in, imgW, imgN, imgT, imgEq, xGrid, yGrid, imgC, imgP, imgPf;
    int i, j;
    int ncs, nrs;

    float *psharp, *pblur;
    float xmin, xmax, ymin, ymax;

    float fcx, fcy;

    float *cblur, *csharp;
    float ps;

    float *A, *b, *x;
    int ncol, nrow;

    ImageFloat imgMask;

    ThinPlate tp, tpI;

    /* convert input image to ImageFloat */
    in = new_imageFloat (nx, ny);
    memcpy (in->val, img, nx * ny * sizeof (float));

    /* Convert the input random pattern image to a sharp pattern image
     of UP_RES x UP_RES larger size by replacing each pixel by a block of
     UP_RES x UP_RES pixels with the same gray value. Also normalize
     the sharp pattern image to be a FloatImage in [0,1]*/
    imgP =  pattern_to_pattern_image(pattern, pat_nx, pat_ny);


    /*---------Detecting the pattern---------*/
    printf ("Detecting the pattern...\n");
    pblur = detect_pattern (in);


    /*Check if the detected pattern is required by the user*/
    if(detected_ppm)
    {
        ImageFloat inR_detected, inG_detected, inB_detected;
        inR_detected =  draw_detected_corners_image_maxval(in, pblur);
        inG_detected =  draw_detected_corners_image_minval(in, pblur);
        inB_detected =  draw_detected_corners_image_minval(in, pblur);

        write_ppm_normalize_float(detected_ppm,
                                  inR_detected->val,
                                  inG_detected->val,
                                  inB_detected->val,
                                  inR_detected->ncol, inR_detected->nrow);
        free_imageFloat(inR_detected);
        free_imageFloat(inG_detected);
        free_imageFloat(inB_detected);

    }

    /*---Extract a sub image with the target
     *  and update the checkpoints locations relative to the extracted image
     */
    printf ("Extracting the pattern...\n");
    imgT = extract_pattern_region (in, pblur, 12);


#ifdef SAVE_INTERMEDIATE_IMGS
    write_pgm_normalize_float("pattern_region_image.pgm",
                              imgT->val,imgT->ncol, imgT->nrow);
#endif

    /*----Geometric Distortion - Thin Plates----*/
    printf ("Calculating thin-plates distortion...\n");

    /* Compute checkerboard X points positions in the analytic pattern */
    psharp = pattern_Xpoints ();

    tp = calculate_thinPlate (pblur, psharp, 12, 10);/*lambda = 10 */
    tpI = calculate_thinPlate (psharp, pblur, 12, 10);/*lambda = 10 */

    /*---Image Ilumination Normalization B&W---*/
    printf ("Normalizing image illumination...\n");
    imgN = image_normalization (imgT, tpI);
#ifdef SAVE_INTERMEDIATE_IMGS
    write_pgm_normalize_float("illumination_normalized_image.pgm",
                              imgN->val, imgN->ncol, imgN->nrow);
#endif

    /*---CRF Estimation and Correction---*/
    printf ("Estimating and correcting CRF...\n");
    imgEq = crf_correction (imgN, tpI);
#ifdef SAVE_INTERMEDIATE_IMGS
    write_pgm_normalize_float("crf_corrected_image.pgm",
                              imgEq->val, imgEq->ncol, imgEq->nrow);
#endif

    /*---Extract Noise Region from the observed image*/
    printf ("Extracting noise region...\n");
    imgW =
    extract_noise_region (imgEq, tpI, &xmin, &xmax, &ymin, &ymax, &imgMask);
#ifdef SAVE_INTERMEDIATE_IMGS
    write_pgm_normalize_float("noise_region_image.pgm",
                              imgW->val, imgW->ncol, imgW->nrow);
    write_pgm_normalize_float("mask_image.pgm",
                              imgMask->val, imgMask->ncol, imgMask->nrow);
#endif

    /*--- Pattern Rasterization --- */
    /*--------> Cut the spectrum- */
    /*
     * int m = up_res*512;
     * int n = up_res*512;
     * double fcx = q*s/(2*n);
     * double fcy = p*s/(2*m);
     */

    /*
     fcx = imgW->ncol*s/(2*UP_RES*512);
     fcy = imgW->nrow*s/(2*UP_RES*512);
     */
    printf ("Filtering and interpolating sharp pattern image...\n");

    fcx = roundfi (imgW->ncol * s);
    fcy = roundfi (imgW->nrow * s);

    imgPf = lpf_image_dct (imgP, (int) fcx, (int) fcy);


    /* Pattern sx Interpolation */
    ps = 1 / ((float) s);
    ncs = (xmax - xmin) * s + 1;
    nrs = (ymax - ymin) * s + 1;

    cblur = (float *) malloc (nrs * ncs * 2 * sizeof (float));
    csharp = (float *) malloc (nrs * ncs * 2 * sizeof (float));

    /*Generating the sx-sampling grid */
    for (i = 0; i < nrs; i++)
        for (j = 0; j < ncs; j++)
        {
            cblur[2 * i * ncs + 2 * j] = xmin + j * ps;
            cblur[2 * i * ncs + 2 * j + 1] = ymin + i * ps;

        }

    /*Applying TP to the sx-sampling grid and interpolate the Sharp pattern */
    evaluate_thinPlate (tp, cblur, csharp, nrs * ncs);

    xGrid = new_imageFloat (ncs, nrs);
    yGrid = new_imageFloat (ncs, nrs);

    /*The Sharp pattern image only contains the noise part
     * so i need to translate the x,y Grid. Noise regions
     * starts in (PATTERN_BLOCK_SIZE,PATTERN_BLOCK_SIZE)*UP_RES
     */
    for (i = 0; i < nrs; i++)
        for (j = 0; j < ncs; j++)
        {
            xGrid->val[i * ncs + j] =
            csharp[2 * i * ncs + 2 * j] - PATTERN_BLOCK_SIZE * UP_RES;
            yGrid->val[i * ncs + j] =
            csharp[2 * i * ncs + 2 * j + 1] - PATTERN_BLOCK_SIZE * UP_RES;

        }

    /*In version 1.0 the parameter was wrongly set to 0.5 instead of -0.5*/
    imgC = bicubic (xGrid, yGrid, imgPf, -0.5);
#ifdef SAVE_INTERMEDIATE_IMGS
    write_pgm_normalize_float("pattern_interpolated_image.pgm",
                              imgC->val, imgC->ncol, imgC->nrow);
#endif

    /*Generating A and b for Ax = b system */
    printf ("Generating A,b for Ax = b linear system...\n");
    make_Ab (imgC, imgW, imgMask, psf_nrow, psf_ncol, s, &A, &b, &ncol, &nrow);


    /*There are three different ways of solving
     * x / Ax = b.
     * i)   least squares
     * ii)  least squares and then projection (x>=th)
     * iii) non-negative least squares (x>=0)
     */

    printf ("Estimating the PSF: solving Ax = b...\n");

    if(solver==0)
    {
        x = solve_lsd(A, b, nrow, ncol);
    }
    else if (solver==1)
    {
        x = solve_lsd_th (A, b, nrow, ncol, 0.001);
    }
    else if (solver ==2)
    {
        x = solve_nnlsd(A, b, nrow, ncol);
    }
    else {
        error("Solver not recognized. Solver must be 0,1,2");
    }

    printf ("Cleaning the house...\n");
    free_imageFloat (imgW);
    free_imageFloat (imgN);
    free_imageFloat (imgT);
    free_imageFloat (imgEq);
    free_imageFloat (xGrid);
    free_imageFloat (yGrid);
    free_imageFloat (imgC);
    free_imageFloat (imgP);
    free_imageFloat (imgPf);
    free_imageFloat (in);

    free_thinPlate (tp);
    free_thinPlate (tpI);

    free ((void *) cblur);
    free ((void *) csharp);
    free ((void *) psharp);
    free ((void *) pblur);
    free ((void *) A);
    free((void*) b);


    return x;
}
