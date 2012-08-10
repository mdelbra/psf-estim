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
 * @file detect_pattern.c
 * @brief  Module code for detecting the pattern in a digital image
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "image.h"
#include "lsd.h"
#include "detect_pattern.h"

/** Number of segments that will be detected in the pattern image  */
#define NUM_SEG 13

/** Relative Tolerance to consider a segment is detected */
#define TOL 0.1

/** Big Number */
#define BIG_NUMBER 1000000

/** Maximum LSD scale */
#define MAX_SCALE_LSD 4.0

/** Initial LSD scale */
#define INITIAL_SCALE_LSD  1.85

static void error(char *msg)
{
    fprintf(stderr, "Detect_Pattern Error: %s\n", msg);
    exit(EXIT_FAILURE);
}



/* homemade round function: float to int*/
int roundfi(float x)
{
    if ((x <= INT_MIN-0.5) || (x >= INT_MAX+0.5))
        error("roundfi() Float to int conversion out of range");
    if (x >= 0)
        return (int) (x+0.5);
    return (int) (x-0.5);
}


/**
 * @brief Coarse detection of the pattern by using segments detected from LSD
 * @param in - input image float
 * @return Array of floats containing the OPQR 2D point positions of the
 *         pattern
 */
static float *detect_pattern_coarse(ImageFloat in)
{
    image_double in_lsd;
    ntuple_list seg;
    float min_val, max_val;
    double *seg_length, *seg_midpoint_x, *seg_midpoint_y;
    char is_in_distance;
    double xO, xP, xQ, xR, yO, yQ, yP, yR, l;
    int *has_seg, *has_all_seg, point;

    /* To keep track of the error from the optimal segment position
     * just to choose the closest segment to the optimal location
     */
    double seg_error[NUM_SEG];

    double xi1, xi2, yi1, yi2, xj1, xj2, yj1, yj2, xjN1, xjN2, yjN1, yjN2;
    double xoMp, xpMq, xqMr, xrMo, yoMp, ypMq, yqMr, yrMo;
    double xC, yC;

    /* Structure of the Pattern  Seg0,Seg1,...,Seg10,Seg_Orient0,Seg_Orient1 */
    const double x1S[] = { 3, 6, 7, 7, 7, 6, 3, 0, -2, -2, -2, 1, 8 };
    const double y1S[] = { 0, 0, -1, -4, -7, -9, -9, -9, -7, -4, -1, 1, 0 };
    const double x2S[] = { 3, 6, 8, 8, 8, 6, 3, 0, -1, -1, -1, 1, 7 };
    const double y2S[] = { 1, 1, -1, -4, -7, -8, -8, -8, -7, -4, -1, 0, 0 };


    int i = 0, j = 0, k = 0, c = 0;

    double errorc;
    double actual_scale;
    char ready;

    float *opqr = (float *) malloc(8 * sizeof(float));


    /* Execute LSD */
    /* Convert between images types and renormalize the image to [0,255]*/
    min_val = BIG_NUMBER;
    max_val = 0;
    for (i=0; i< in->ncol *in->nrow;i++)
    {
        if(in->val[i] < min_val) min_val = in->val[i];
        if(in->val[i] > max_val) max_val = in->val[i];
    }

    in_lsd = new_image_double((unsigned int) in->ncol,
                               (unsigned int) in->nrow);
    for (i = 0; i < in->ncol * in->nrow; i++)
        in_lsd->data[i] = (double) 255/(max_val-min_val)*(in->val[i]-min_val);

    /* We do a LOOP from INITIAL_SCALE_LSD to MAX_SCALE_LSD */
    actual_scale = INITIAL_SCALE_LSD;
    ready = 0;

    while (actual_scale < MAX_SCALE_LSD && !ready)
    {
        printf(" -->LSD scale =%f\n",actual_scale);
        seg = lsd_scale(in_lsd, actual_scale);

        /* allocate the array or exit with an error */
        if ((seg_length = (double *) malloc(seg->size * sizeof(double)))
                  == NULL
            || (seg_midpoint_x = (double *) malloc(seg->size * sizeof(double)))
                  == NULL
            || (seg_midpoint_y = (double *) malloc(seg->size * sizeof(double)))
                  == NULL)
        {
            error("PSF_ESTIM - Unable to allocate double array space");
            exit(EXIT_FAILURE);
        }

        /*
         The i component, of the n-tuple number j, of an n-tuple list 'ntl'
         is accessed with:
         */
        for (i = 0; i < (int) seg->size; i++)
        {

            /* segment length */
            seg_length[i] = dist_l2(seg->values[i * seg->dim],
                  seg->values[i * seg->dim + 1],
                  seg->values[i * seg->dim + 2],
                  seg->values[i * seg->dim + 3]);

            /* segment midpoint */
            seg_midpoint_x[i] =
                0.5 * (seg->values[i * seg->dim]
                + seg->values[i * seg->dim + 2]);

            seg_midpoint_y[i] =
                0.5 * (seg->values[i * seg->dim + 1]
                + seg->values[i * seg->dim + 3]);
        }

        /* Accessing to segment j=0...12 associated to segment i
         * has_seg[NUM_SEG*i+j], initialization default to 0
         */

        if ((has_seg =
             (int *) malloc(NUM_SEG * seg->size * sizeof(int))) == NULL
            || (has_all_seg =
                (int *) malloc(seg->size * sizeof(int))) == NULL)
        {
            error("PSF_ESTIM - Unable to allocate double array space");
            exit(EXIT_FAILURE);
        }

        /* has_seg[], Initialization default to -1 */
        for (i = 0; i < NUM_SEG * (int) seg->size; i++)
        {
            has_seg[i] = -1;
        }

        /* has_all_seg[], Initialization default to 0 */
        /* First pass */
        for (i = 0; i < (int) seg->size; i++)
        {
            xi1 = seg->values[i * seg->dim];
            xi2 = seg->values[i * seg->dim + 2];
            yi1 = seg->values[i * seg->dim + 1];
            yi2 = seg->values[i * seg->dim + 3];

            /* Reinitialize the error track */
            for (j = 0; j < NUM_SEG;  j++)
            {
                seg_error[j] = BIG_NUMBER;
            }

            for (j = 0; j < (int) seg->size; j++)
            {
                xj1 = seg->values[j * seg->dim];
                xj2 = seg->values[j * seg->dim + 2];
                yj1 = seg->values[j * seg->dim + 1];
                yj2 = seg->values[j * seg->dim + 3];

                /* Convert the (x,y) coordinates to a new Coordinate System
                 * (xN, yN) having:
                 * (xi1,yi1) at (0,0)
                 * (xi2,yi2) at (0,1)
                 *  The vectors x1s,y1s,x2s,y2s are given within this
                 *  new (xN,yN) coordinate system
                 */
                l = seg_length[i];
                xjN1 =
                    1 / (l * l) * ((yi2 - yi1) * (xj1 - xi1)
                                   - (xi2 - xi1) * (yj1 - yi1));
                yjN1 =
                    1 / (l * l) * ((xi2 - xi1) * (xj1 - xi1)
                                   + (yi2 - yi1) * (yj1 - yi1));
                xjN2 =
                    1 / (l * l) * ((yi2 - yi1) * (xj2 - xi1)
                                   - (xi2 - xi1) * (yj2 - yi1));
                yjN2 =
                    1 / (l * l) * ((xi2 - xi1) * (xj2 - xi1)
                                   + (yi2 - yi1) * (yj2 - yi1));

                for (c = 0; c < NUM_SEG; c++)
                {
                    is_in_distance =
                        (fabs(xjN1 - x1S[c]) < TOL * (2 + fabs(x1S[c])))
                        && (fabs(yjN1 - y1S[c]) < TOL * (2 + fabs(y1S[c])))
                        && (fabs(xjN2 - x2S[c]) < TOL * (2 + fabs(x2S[c])))
                        && (fabs(yjN2 - y2S[c]) < TOL * (2 + fabs(y2S[c])));
                    if (is_in_distance)
                    {
                        /* Need to check that there isn't a previous segment
                         * closer to the optimal location and already marked
                         * as good (errorc). I just keep the segment with
                         * minimum total error l1.
                         */
                        errorc = fabs(xjN1 - x1S[c]) + fabs(yjN1 - y1S[c])
                            + fabs(xjN2 - x2S[c]) + fabs(yjN2 - y2S[c]);
                        if (errorc < seg_error[c])
                        {
                            has_seg[i * NUM_SEG + c] = j;
                            seg_error[c] = errorc;
                        }
                    }
                }
            }


            /*has_all_seg[i] will be one if all segments are present */
            has_all_seg[i] = 1;
            for (j=0;j< NUM_SEG;j++)
            {
                has_all_seg[i] =
                     has_all_seg[i] && (has_seg[i * NUM_SEG + j] >= 0);
            }

            if (has_all_seg[i])
            {
                point = i;
                k++;
            }
        }
        ready = (k==1);
        actual_scale *= 1.15;
    }

    if (k > 1)
    {
        printf("More than one pattern was detected.");
        printf("\nCrop the image surounding the desired pattern and re-run.");
        exit(EXIT_SEVERAL_PATTERNS_DETECTED);
    }

    else if(k<1)
    {
        printf("No pattern was detected. Use another image");
        exit(EXIT_NO_PATTERN_DETECTED);
    }


    /*Calculate C - center, u = unit_length, theta = angle */
    /* 1/3*(DET + 0 + 1) = oMp */
    xoMp = 0.33333 * (seg_midpoint_x[point]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 0]]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 1]]);

    yoMp = 0.33333 * (seg_midpoint_y[point]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 0]]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 1]]);

    /* 1/3*(2 + 3 + 4) = pMq */
    xpMq = 0.33333 * (seg_midpoint_x[has_seg[point * NUM_SEG + 2]]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 3]]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 4]]);

    ypMq = 0.33333 * (seg_midpoint_y[has_seg[point * NUM_SEG + 2]]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 3]]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 4]]);

    /* 1/3*(5 + 6 + 7) = qMr */
    xqMr = 0.33333 * (seg_midpoint_x[has_seg[point * NUM_SEG + 5]]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 6]]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 7]]);

    yqMr = 0.33333 * (seg_midpoint_y[has_seg[point * NUM_SEG + 5]]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 6]]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 7]]);

    /* 1/3*(8 + 9 + 10) = rMo */
    xrMo = 0.33333 * (seg_midpoint_x[has_seg[point * NUM_SEG + 8]]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 9]]
                      + seg_midpoint_x[has_seg[point * NUM_SEG + 10]]);

    yrMo = 0.33333 * (seg_midpoint_y[has_seg[point * NUM_SEG + 8]]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 9]]
                      + seg_midpoint_y[has_seg[point * NUM_SEG + 10]]);

    /*Center */
    xC = 0.25 * (xoMp + xpMq + xqMr + xrMo);
    yC = 0.25 * (yoMp + ypMq + yqMr + yrMo);

    /*O = C + CoMr + CoMp */
    xO = xC + (xrMo - xC) + (xoMp - xC);
    yO = yC + (yrMo - yC) + (yoMp - yC);

    /*P = C + CpMq + CoMp */
    xP = xC + (xpMq - xC) + (xoMp - xC);
    yP = yC + (ypMq - yC) + (yoMp - yC);

    /*Q = C + CqMr + CpMq */
    xQ = xC + (xqMr - xC) + (xpMq - xC);
    yQ = yC + (yqMr - yC) + (ypMq - yC);

    /*R = C + CrMo + CqMr */
    xR = xC + (xrMo - xC) + (xqMr - xC);
    yR = yC + (yrMo - yC) + (yqMr - yC);

    /*Array of OPQR coordinates*/
    opqr[0] = (float) xO;
    opqr[1] = (float) yO;
    opqr[2] = (float) xP;
    opqr[3] = (float) yP;
    opqr[4] = (float) xQ;
    opqr[5] = (float) yQ;
    opqr[6] = (float) xR;
    opqr[7] = (float) yR;

    /* free memory */
    free_image_double(in_lsd);

    free_ntuple_list(seg);
    free(seg_length);
    free(seg_midpoint_y);
    free(seg_midpoint_x);
    free(has_seg);
    free(has_all_seg);

    return opqr;
}




/**
 * @brief Detection of a X corner in the imput image
 * @param in - input image float
 * @param ptx - approximate x coord. of the X corner; (output) refined position
 * @param pty - approximate y coord. of the X corner; (output) refined position
 * @return int - 0 if no error
 */
static int detect_xcorner(ImageFloat in, float *ptx, float *pty)
{

    ImageFloat mask, src_buff, gx_buff, gy_buff;
    float coeff;
    int i, j, k, y, x;
    float cx = *ptx;
    float cy = *pty;
    float c2x, c2y;

    int max_iters = 400;

    float eps = 0.00001;
    int wsize = 3;
    float a11, a12, a22, p, q, d;
    int iter = 0;
    float err;
    float py, px;
    float tgx, tgy, gxx, gxy, gyy, m;

    float *mask1D;


    /*mask1D = new_vector(2*wsize+1); */

    mask1D = (float *) malloc((2 * wsize + 1) * sizeof(float));

    mask = new_imageFloat(2 * wsize + 1, 2 * wsize + 1);

    coeff = 1. / (mask->ncol * mask->nrow);

    /* calculate mask */
    for (i = -wsize, k = 0; i <= wsize; i++, k++)
    {
        mask1D[k] = exp(-i * i * coeff);
    }

    for (i = 0; i < (int) mask->nrow; i++)
    {
        for (j = 0; j < (int) mask->ncol; j++)
        {
            mask->val[i * mask->nrow + j] = mask1D[j] * mask1D[i];
        }
    }


    do {
        src_buff = extract_subpx_window(in, wsize, cx, cy);
        gx_buff = gradx(src_buff);
        gy_buff = grady(src_buff);
        a11 = a12 = a22 = p = q = 0;

        /* process gradient */
        for (y = -wsize, k = 0; y <= wsize; y++)
        {
            py = cy + (float) y;

            for (x = -wsize; x <= wsize; x++, k++)
            {
                m = mask->val[k];
                tgx = gx_buff->val[k];
                tgy = gy_buff->val[k];
                gxx = tgx * tgx * m;
                gxy = tgx * tgy * m;
                gyy = tgy * tgy * m;
                px = cx + (float) x;

                a11 += gxx;
                a12 += gxy;
                a22 += gyy;

                p += gxx * px + gxy * py;
                q += gxy * px + gyy * py;
            }
        }

        d = a11 * a22 - a12 * a12;
        c2x = 1 / d * (a22 * p - a12 * q);
        c2y = 1 / d * (-a12 * p + a11 * q);

        err = dist_l2(cx, cy, c2x, c2y);
        cx = c2x;
        cy = c2y;

        free_imageFloat(src_buff);
        free_imageFloat(gx_buff);
        free_imageFloat(gy_buff);
    } while (++iter < max_iters && err > eps);


    *ptx = cx;
    *pty = cy;

    free_imageFloat(mask);
    free((void *) mask1D);

    return 0;
}


/**
 * @brief Precise detection of the pattern by using the Coarse estimation
 * @param in - input image float
 * @param opqr - Array of floats containing the OPQR 2D point positions of
 *        the pattern
 * @return Array of 12 points where the X marks are subpixecally located.
 */
static float *detect_pattern_fine(ImageFloat in, float *opqr)
{
    /* there should be 12 X-checkerboard corners:
     *   O O1 O2 P P1 P2 Q Q1 Q2 R R1 R2
     *0  1  2 3 4  5  6 7   8 9 10 11
     *
     *NOTE: O - is the X corner neighbor of two black squares
     */
    float *p = (float *) malloc(12 * 2 * sizeof(float));

    /* O - initial guess */
    p[0] = opqr[0];
    p[1] = opqr[1];

    /* P - initial guess */
    p[6] = opqr[2];
    p[7] = opqr[3];

    /* Q - initial guess */
    p[12] = opqr[4];
    p[13] = opqr[5];

    /* R - initial guess */
    p[18] = opqr[6];
    p[19] = opqr[7];


    /* first detect opqr at subpixel precision */
    detect_xcorner(in, &p[0], &p[1]);
    detect_xcorner(in, &p[6], &p[7]);
    detect_xcorner(in, &p[12], &p[13]);
    detect_xcorner(in, &p[18], &p[19]);

    /*linear interpolation as initial position guess
     * of the secondary points
     */

    /* O1 - initial guess and subpixel detection */
    p[2] = 0.6666 * p[0] + 0.3333 * p[6];
    p[3] = 0.6666 * p[1] + 0.3333 * p[7];
    detect_xcorner(in, &p[2], &p[3]);

    /* O2 - initial guess and subpixel detection */
    p[4] = 0.3333 * p[0] + 0.6666 * p[6];
    p[5] = 0.3333 * p[1] + 0.6666 * p[7];
    detect_xcorner(in, &p[4], &p[5]);


    /* P1 - initial guess and subpixel detection */
    p[8] = 0.6666 * p[6] + 0.3333 * p[12];
    p[9] = 0.6666 * p[7] + 0.3333 * p[13];
    detect_xcorner(in, &p[8], &p[9]);

    /* P2 - initial guess and subpixel detection */
    p[10] = 0.3333 * p[6] + 0.6666 * p[12];
    p[11] = 0.3333 * p[7] + 0.6666 * p[13];
    detect_xcorner(in, &p[10], &p[11]);


    /* Q1 - initial guess and subpixel detection */
    p[14] = 0.6666 * p[12] + 0.3333 * p[18];
    p[15] = 0.6666 * p[13] + 0.3333 * p[19];
    detect_xcorner(in, &p[14], &p[15]);

    /* Q2 - initial guess and subpixel detection */
    p[16] = 0.3333 * p[12] + 0.6666 * p[18];
    p[17] = 0.3333 * p[13] + 0.6666 * p[19];
    detect_xcorner(in, &p[16], &p[17]);


    /* R1 - initial guess and subpixel detection */
    p[20] = 0.6666 * p[18] + 0.3333 * p[0];
    p[21] = 0.6666 * p[19] + 0.3333 * p[1];
    detect_xcorner(in, &p[20], &p[21]);

    /* R2 - initial guess and subpixel detection */
    p[22] = 0.3333 * p[18] + 0.6666 * p[0];
    p[23] = 0.3333 * p[19] + 0.6666 * p[1];
    detect_xcorner(in, &p[22], &p[23]);



    return p;

}

/**
 * @brief  Detection of the pattern by using the Coarse & Precise estimations
 *    There is a fisrt Coarse estimation by using LSD and then a second pass
 *    in order to refine the position of the X-checkerboard marks presented in
 *     the pattern
 * @param in - input image float
 * @return Array of 12 points where the X marks are subpixecally located.
 */
float *detect_pattern(ImageFloat in)
{

    float *opqr, *checkpoints;

    /*opqr pattern location */

    /*X-Checkpoints */
    /*O O1 O2 P P1 P2 Q Q1 Q2 R R1 R2 pattern location */

    /* coarse detection of the pattern by using LSD */
    opqr = detect_pattern_coarse(in);

    /* precise detection of the pattern by using the checkerboard marks
     and Bouguet-OpenCV X-detector */
    checkpoints = detect_pattern_fine(in, opqr);

    /*Clean opqr */
    free((void *) opqr);
    
    return checkpoints;

}


/**
 * @brief  Gives the positions of the X marks in the analityc pattern
 * @return Array of 12 points where the X marks are subpixecally located.
 */
float *pattern_Xpoints(void)
{
    /*This function returns the locations - in the analytic pattern -
     * where the X checkerboard corners are.
     */

    float *p = (float *) malloc(12 * 2 * sizeof(float));

    /* 12 points O..P..Q..R.. coord x (odd) and y (even) */
    float up_res = UP_RES;

    int up_res_pattern_block_size = up_res * PATTERN_BLOCK_SIZE;
    int up_res_pattern_block_size_2 = up_res_pattern_block_size/2;


    /*From Point O in counterclockwise order*/
    /*O*/
    p[0] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[1] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    /*O1 */
    p[2] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[3] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 6;

    /*O2 */
    p[4] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[5] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 3;

    /*P*/
    p[6] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[7] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*P1 */
    p[8] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 6;
    p[9] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*P2 */
    p[10] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 3;
    p[11] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*Q*/
    p[12] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[13] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*Q1 */
    p[14] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[15] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 3;

    /*Q2 */
    p[16] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[17] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 6;

    /*R*/
    p[18] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[19] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    /*R1 */
    p[20] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 3;
    p[21] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    /*R2 */
    p[22] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 6;
    p[23] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    return p;


}

/**
 * @brief  Gives the positions of the black squares centers in the analytic
 *        pattern
 * @param p Array of 12 points where the black squares centers are
 *         subpixecally located on output
 */
void pattern_blackSquares(float *p)
{
    /* This function returns the locations - in the analytic pattern -
     * where the center of the black squares are placed.
     */


    /* 12 points  coord x (odd) and y (even) */
    float up_res = UP_RES;

    int up_res_pattern_block_size = up_res * PATTERN_BLOCK_SIZE;
    int up_res_pattern_block_size_2 = up_res_pattern_block_size/2;


    /*From Point O in counterclockwise order */
    /*O*/
    p[0] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[1] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 8;

    /*O1 */
    p[2] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[3] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 5;

    /*O2 */
    p[4] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[5] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 2;

    /*P*/
    p[6] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 8;
    p[7] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*P1 */
    p[8] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 5;
    p[9] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*P2 */
    p[10] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 2;
    p[11] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*Q*/
    p[12] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[13] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 2;

    /*Q1 */
    p[14] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[15] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 5;

    /*Q2 */
    p[16] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[17] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 8;

    /*R*/
    p[18] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 2;
    p[19] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    /*R1 */
    p[20] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 5;
    p[21] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    /*R2 */
    p[22] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 8;
    p[23] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    return;

}

/**
 * @brief  Gives the positions of the black squares centers in the analytic
 *        pattern
 * @param p Array of 12 points where the white squares centers are
 *        subpixecally located on output
 */
void pattern_whiteSquares(float *p)
{
    /*This function returns the locations - in the analytic pattern -
     * where the center of the white squares are placed.
     */


    /* 12 points  coord x (odd) and y (even) */
    float up_res = UP_RES;

    int up_res_pattern_block_size = up_res * PATTERN_BLOCK_SIZE;
    int up_res_pattern_block_size_2 = up_res_pattern_block_size/2;


    /*From Point O in counterclockwise order */
    /*O*/
    p[0] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[1] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 7;

    /*O1 */
    p[2] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[3] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 4;

    /*O2 */
    p[4] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;
    p[5] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 1;

    /*P*/
    p[6] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 7;
    p[7] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*P1 */
    p[8] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 4;
    p[9] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*P2 */
    p[10] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 1;
    p[11] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;

    /*Q*/
    p[12] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[13] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 1;

    /*Q1 */
    p[14] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[15] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 4;

    /*Q2 */
    p[16] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 0;
    p[17] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 7;

    /*R*/
    p[18] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 1;
    p[19] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    /*R1 */
    p[20] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 4;
    p[21] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    /*R2 */
    p[22] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 7;
    p[23] = up_res_pattern_block_size_2 - 0.5 + up_res_pattern_block_size * 9;

    return;

}

/**
 * @brief  Gives the positions of the the center point inside the noise region
 * @param p Array of 1 point where the point is located
 */
void pattern_center(float *p)
{
    /*This function returns the location - in the analytic pattern -
     * of the pattern center point.
     */

    /* 1 points  coord x (odd) and y (even) */
    float up_res = UP_RES;

    p[0] = PATTERN_BLOCK_SIZE * up_res * 5 - 0.5;
    p[1] = PATTERN_BLOCK_SIZE * up_res * 5 - 0.5;

    return;

}

/**
 * @brief  Gives the positions of the the topmost point inside the noise region
 *          part at the horizontal center
 * @param p Array of 1 point where the point is located
 */
void pattern_top_center(float *p)
{

    /* 1 points  coord x (odd) and y (even) */
    float up_res = UP_RES;

    p[0] = PATTERN_BLOCK_SIZE * up_res * 5 - 0.5;
    p[1] = PATTERN_BLOCK_SIZE * up_res * 1 - 0.5;

    return;

}


/**
 * @brief  Draw a X at position 'x','y' of length 'w' on image 'in'
 * @param x horizontal coordiante (integer)
 * @param y vertical coordiante (integer)
 * @param w length of the X-mark (integer)
 */
static void draw_x(int x, int y, int w, float val, ImageFloat in)
{
    int i;
    for(i=-w;i<w;i++)
    {
        in->val[x + i + (y+i)*in->ncol] = val;
        in->val[x - i -1 + (y+i)*in->ncol] = val;

    }

    return;
}

/**
 * @brief  Draw the detected corners with intensity equal to  the max value
 *         in the input image
 * @param in input image float
 * @param checkerboard arrat containing the 12 locations of the X-checkerboards
 * @return float image with X's where the points are located
 */
ImageFloat draw_detected_corners_image_maxval(ImageFloat in,
                                              float *checkerboard)
{
    /*Checkpoints */
    /*O O1 O2 P P1 P2 Q Q1 Q2 R R1 R2 pattern location */
    int length = 4;
    int  i;
    float max_val=0;
    ImageFloat in_detected;

    in_detected = new_imageFloat(in->ncol,in->nrow);

    for(i=0;i<in->nrow*in->ncol;i++)
    {
        in_detected->val[i] = in->val[i];
        if(in->val[i]>max_val) max_val = in->val[i];
    }

    for(i=0; i<12;i++)
        /*0.5 is added to draw the segment in the center of the pixel*/
        draw_x(roundfi(checkerboard[2*i]+0.5),
               roundfi(checkerboard[2*i+1]+0.5),
               length, max_val, in_detected);

    return in_detected;
}

/**
 * @brief  Draw the detected corners with intensity equal to  the min value
 *         in the input image
 * @param in input image float
 * @param checkerboard arrat containing the 12 locations of the X-checkerboards
 * @return float image with X's where the points are located
 */
ImageFloat draw_detected_corners_image_minval(ImageFloat in,
                                              float *checkerboard)
{
    /*Checkpoints */
    /*O O1 O2 P P1 P2 Q Q1 Q2 R R1 R2 pattern location */
    int length = 4;
    int  i;
    float min_val = BIG_NUMBER;
    ImageFloat in_detected;

    in_detected = new_imageFloat(in->ncol,in->nrow);

    for(i=0;i<in->nrow*in->ncol;i++)
    {
        in_detected->val[i] = in->val[i];
        if(in->val[i]<min_val) min_val = in->val[i];
    }

    for(i=0; i<12;i++)
        /*0.5 is added to draw the segment in the center of the pixel*/
        draw_x(roundfi(checkerboard[2*i]+0.5),
                       roundfi(checkerboard[2*i+1]+0.5),
                       length, min_val, in_detected);

    return in_detected;
}


/**
 * @brief Convert the input random pattern image to a sharp pattern image
          of UP_RES x UP_RES larger size by replacing each pixel by a block of
          UP_RES x UP_RES pixels with the same gray value. Also normalize
          the sharp pattern image to be a FloatImage in [0,1]
 * @param pattern pattern float
 * @param pat_nx horizontal size of the pattern input
 * @param pat_nx vertical size of the pattern input
 * @return ImageFloat with pattern rasterized at UP_RES resolution
 */
ImageFloat pattern_to_pattern_image(float* pattern, int pat_nx, int pat_ny)
{
    int maxval, i, j, k, l;
    float pixval;
    ImageFloat imgP;

    maxval = 0;

    for(i=0;i<pat_nx*pat_ny;i++)
        if(pattern[i]>maxval) maxval = pattern[i];

    imgP = new_imageFloat ((int) pat_nx * UP_RES, (int) pat_ny*UP_RES);
    for(i=0; i < pat_ny;i++)
        for (j=0;j< pat_nx;j++)
        {
            pixval = pattern[j + i*pat_nx]/maxval;
            for(k=0;k< (int) UP_RES;k++)
                for(l=0;l< (int) UP_RES;l++)
                    imgP->val[j* (int)UP_RES
                              + l + imgP->ncol*(i* (int)UP_RES + k)]
                    = pixval;

        }


    return imgP;
}
