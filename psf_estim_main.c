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
 * @file psf_estim_main.c
 * @brief main for psf estimation algorithm execution.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "io_pgm.h"
#include "image.h"
#include "psf_estim.h"


/** @brief Struct of program parameters */
typedef struct
{
    int s;
    int psf_nx;
    int psf_ny;
    char *pattern_pgm;
    char *out_pgm;
    char *input;
    char *output;
    char *detected_ppm;
    int solver;
} program_argums;


/** @brief Error/Exit print a message and exit.
 *  @param msg
 */
static void error(char *msg)
{
    fprintf(stderr, "PSF_ESTIM Error: %s\n", msg);
    exit(EXIT_FAILURE);
}


static void usage(const char* name)
{
    printf("Point Spread Function Estimation from a Random Target\n");
    printf("Copyright M.Delbracio, P.Muse, A.Almansa. ");
    printf("Version 1.3 - 26 December 2011\n\n");
    printf("Usage: %s [options] <input file> <output file>\n\n"
           "Only  PGM 16/8 bits images are supported.\n\n",name);
    printf("Options:\n");
    printf("  -s <number>  The super-resolution factor, positive integer");
    printf("  (default 4)\n");
    printf("  -k <number>  PSF support size (default 4s + 1)\n");
    printf("  -o <file>    Estimated PSF save to a 8bits PGM image \n");
    printf("  -p <file>    Noise Pattern PGM image \n");
    printf("  -d <file>    Save the input image with detected pattern marked");
    printf("  (24bits PPM)\n");
    printf("  -t <0,1,2>   Least Squares solver:\n");
    printf("                    0 - Least Squares\n");
    printf("                    1 - Least Squares + thresholding (default)\n");
    printf("                    2 - Non-negative least Squares (slowest)\n");
}

static void parse_arguments(program_argums *param, int argc, char *argv[])
{
    char *OptionString;
    char OptionChar;
    int i;


    if(argc < 2)
    {
        usage(argv[0]);
        exit(EXIT_SUCCESS);
    }


    /* loop to read parameters*/
    for(i = 1; i < argc;)
    {
        if(argv[i] && argv[i][0] == '-')
        {
            if((OptionChar = argv[i][1]) == 0)
            {
                error("Invalid parameter format.\n");
            }

            if(argv[i][2])
                OptionString = &argv[i][2];
            else if(++i < argc)
                OptionString = argv[i];
            else
            {
                error("Invalid parameter format.\n");
            }

            switch(OptionChar)
            {
                case 's':
                    param->s = atoi(OptionString);
                    if(param->s < 1)
                    {
                        error("Invalid superresolution factor.\n");
                    }
                    break;

                case 'k':
                    param->psf_nx = atoi(OptionString);
                    param->psf_ny = atoi(OptionString);
                    if(param->psf_nx <= 0)
                    {
                        error("Invalid PSF support size.\n");
                    }
                    break;


                case 't':
                    param->solver = atoi(OptionString);
                    if(param->solver != 0 &&
                       param->solver != 1 &&
                       param->solver != 2)
                    {
                        error("t must be 0, 1 or 2.\n");
                    }
                    break;

                case 'o':
                    param->out_pgm = OptionString;
                    break;

                case 'd':
                    param->detected_ppm = OptionString;
                    break;

                case 'p':
                    param->pattern_pgm = OptionString;
                    break;

                case '-':
                    usage(argv[0]);
                    exit(EXIT_FAILURE);

                default:
                    if(isprint(OptionChar))
                    {
                        fprintf(stderr, "Unknown option \"-%c\".\n",
                                OptionChar);
                        exit(EXIT_FAILURE);
                    } else
                        error("Unknown option.\n");
            }

        }
        else
        {
            if(!param->input)
                param->input = argv[i];
            else
                param->output = argv[i];

        }

        i++;
    }

    if(!param->input || !param->output)
    {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }


    /* If parameters weren't set, set deafult parameters*/
    param->s = param->s>0 ? param->s : 4;
    param->psf_nx = param->psf_nx>0 ? param->psf_nx : 4*param->s+1;
    param->psf_ny = param->psf_ny>0 ? param->psf_ny : 4*param->s+1;
    param->solver = param->solver>=0 ? param->solver : 1;
    param->pattern_pgm = param->pattern_pgm ?
                            param->pattern_pgm : "pattern_noise.pgm";



    /*Display selected parameters*/
    printf("\n");
    printf("   Loaded arguments: \n");
    printf("   ----------------- \n");
    printf("         Superresolution        s : %dx\n",param->s);
    printf("         PSF Support            k : %dx%d\n",
           param->psf_nx,param->psf_ny);
    printf("         PSF Output (PGM)       o : %s\n",
           param->out_pgm?param->out_pgm:"(no)");
    printf("         Pattern  (PGM)         p : %s\n",
           param->pattern_pgm?param->pattern_pgm:"(no)");
    printf("         Detected Pattern (PPM) d : %s\n",
           param->detected_ppm?param->detected_ppm:"(no)");
    printf("         LS Solver              t : %d\n",param->solver);
    printf("         PGM input                : %s\n",param->input);
    printf("         PSF Output (TXT)         : %s\n\n",param->output);


}


/**
 * @brief main function call
 */
int main(int argc, char *argv[])
{
    int nx,ny, pat_nx,pat_ny;

    float *in, *img_pattern, *x;

    /*Initialize the structure param->* to -1 or null */
    program_argums param = {-1, -1, -1, NULL, NULL, NULL, NULL, NULL, -1};

    /*Parse command-line arguments*/
    parse_arguments(&param,argc,argv);

    /* read the PGM image into data */
    if (NULL == (in = read_pgm_float(param.input, &nx, &ny))) {
        fprintf(stderr, "the image could not be properly read\n");
        return EXIT_FAILURE;
    }

    /* read the PGM image pattern*/
    if (NULL ==
        (img_pattern = read_pgm_float(param.pattern_pgm, &pat_nx, &pat_ny)))
    {
    fprintf(stderr, "the pattern image could not be properly read\n");
    return EXIT_FAILURE;
    }

    /* Call psf_estimation */
    x = psf_estim(in,  nx,  ny, img_pattern,  pat_nx, pat_ny,
                  param.s, param.psf_nx, param.psf_ny, param.solver,
                  param.detected_ppm);

    /* Write the estimated PSF to a text file */
    write_ascii_matrix(x, param.psf_nx, param.psf_ny, param.output);


    /* Write the estimated PSF to a 8bit-PGM image file if required*/
    if(param.out_pgm)
    {
        /* Normalize and truncate negative values */
        write_pgm_normalize_float(param.out_pgm, x, param.psf_nx,
                                  param.psf_ny);
    }

    /* do the cleaning */
    free(x);
    free(in);
    free(img_pattern);

    return EXIT_SUCCESS;
}
