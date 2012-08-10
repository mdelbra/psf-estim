/*----------------------------------------------------------------------------
 *The followoing functions for reading/writing PGM image files are based on
 *those of LSD - Line Segment Detector on digital images
 *Copyright 2007-2010 rafael grompone von gioi (grompone@gmail.com)
 *The only modification is to allow the reading of 16 bits PGM images
 *(Most significant byte first) according to
 * http://netpbm.sourceforge.net/doc/pgm.html
 */
/*----------------------------------------------------------------------------

 LSD - Line Segment Detector on digital images

 Copyright 2007-2010 rafael grompone von gioi (grompone@gmail.com)

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

/** @file io_pgm.c
 Auxiliary library for reading and writing PGM 8/16 bits images
 @author rafael grompone von gioi (grompone@gmail.com)
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "image.h"
#include "io_pgm.h"


#ifndef FALSE
#define FALSE 0
#endif /* !FALSE */

#ifndef TRUE
#define TRUE 1
#endif /* !TRUE */

/*---------------------------------------------------------------------------*/
/** Fatal error, print a message to standard-error output and exit.
 */
static void error(const char * msg)
{
    fprintf(stderr,"%s\n",msg);
    exit(EXIT_FAILURE);
}

/*---------------------------------------------------------------------------*/
/** Fatal error, print a message and corresponding filename to standard-error
 *  output and exit.
 */

static void errorf(const char * msg, const char *name)
{
    fprintf(stderr,"%s %s.\n",msg,name);
    exit(EXIT_FAILURE);
}


/*---------------------------------------------------------------------------*/
/*------------------------------ PGM image I/O ------------------------------*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/** Skip white characters and comments in a PGM file.
 */
static void skip_whites_and_comments(FILE * f)
{
    int c;
    do
    {
        while(isspace(c=getc(f))); /* skip spaces */
        if(c=='#') /* skip comments */
            while( c!='\n' && c!='\r' && c!=EOF )
                c=getc(f);
    }
    while( c == '#' || isspace(c) );
    if( c != EOF && ungetc(c,f) == EOF )
        error("Error: unable to 'ungetc' while reading PGM file.");
}

/*---------------------------------------------------------------------------*/
/** Read a ASCII number from a PGM file.
 */
static unsigned int get_num(FILE * f)
{
    unsigned int num;
    int c;

    while(isspace(c=getc(f)));
    if(!isdigit(c)) error("Error: corrupted PGM file.");
    num = (unsigned int) (c - '0');
    while( isdigit(c=getc(f)) ) num = 10 * num + c - '0';
    if( c != EOF && ungetc(c,f) == EOF )
        error("Error: unable to 'ungetc' while reading PGM file.");

    return num;
}


/*---------------------------------------------------------------------------*/
/** Read a PGM 8/16 Bits image file into an array of floats.
 */
float *read_pgm_float(const char * fname, int * ncol, int *nrow)
{
    FILE * f;
    int c,bin;
    int depth,i,j;
    float * data;

    /* open file */
    f = fopen(fname,"rb");
    if( f == NULL ) errorf("Error: unable to open input image file ",fname);

    /* read header */
    if( getc(f) != 'P' ) errorf("Error: not a PGM file ",fname);
    if( (c=getc(f)) == '2' ) bin = FALSE;
    else if( c == '5' ) bin = TRUE;
    else errorf("Error: not a PGM file ",fname);
    skip_whites_and_comments(f);
    *ncol = get_num(f);            /* X size */
    skip_whites_and_comments(f);
    *nrow = get_num(f);            /* Y size */
    skip_whites_and_comments(f);
    depth = get_num(f);            /* depth */
    if(depth==0) fprintf(stderr,
                         "Warning: depth=0, probably invalid PGM file\n");
    /* white before data */
    if(!isspace(c=getc(f))) errorf("Error: corrupted PGM file ",fname);

    /* get memory */
    data = (float *) calloc((*ncol) * (*nrow), sizeof(float));
    if (data == NULL)
        error("not enough memory.");

    /* read data */

    /*If the depth is less than 256, it is 1 byte. Otherwise, it is 2 bytes*/
    if(depth<256)
    {
        for(i=0;i<*nrow;i++)
            for(j=0;j<*ncol;j++)
                data[ j + i * (*ncol)] = bin ?
                (float) getc(f) : (float) get_num(f);
    }
    /*16 bits PGM Most significant byte first
     * see http://netpbm.sourceforge.net/doc/pgm.html
     */
    else {
        for(i=0;i<*nrow;i++)
            for(j=0;j<*ncol;j++)
                /*most significant byte first*/
                data[ j + i * (*ncol) ] = bin ? ((float) getc(f)*256) +
                (float)getc(f) : (float) get_num(f);

    }

    /* close file if needed */
    if( f != stdin && fclose(f) == EOF )
        errorf("Error: unable to close file while reading PGM file ",fname);

    return data;
}

/*---------------------------------------------------------------------------*/
/** Write an array of floats into a PGM file.
 */
void write_pgm_float(const char *fname, const float *data, int ncol, int nrow)
{
    FILE * f;
    int i,j;
    int v,max,min;

    /* check min and max values */
    max = min = 0;
    for(i=0; i< nrow; i++)
        for(j=0; j<ncol; j++)
        {
            v = (int) data[ j + i * ncol];
            if( v > max ) max = v;
            if( v < min ) min = v;
        }

    if( min < 0 ) fprintf(stderr,
                  "Warning: negative values in '%s'.\n",
                          fname);
    if( max > 255 ) fprintf(stderr,
                    "Warning: values exceeding 255 in '%s'.\n",
                            fname);

    /* open file */
    if( strcmp(fname,"-") == 0 ) f = stdout;
    else f = fopen(fname,"w");
    if( f == NULL ) errorf("Error: unable to open output image file ",fname);

    /* write header */
    fprintf(f,"P5\n");
    fprintf(f,"%u %u\n",ncol,nrow);
    fprintf(f,"%d\n",255);

    /* write data */
    for(i=0; i<nrow; i++)
        for(j=0; j<ncol; j++)
            fputc((unsigned char) data[j+i*ncol],f);

    /* close file if needed */
    if( f != stdout && fclose(f) == EOF )
        errorf("Error: unable to close file while writing PGM file ",fname);
}


/*---------------------------------------------------------------------------*/
/** Write an array of floats into a 8bit PGM file.
 *    Normalize the maximum value of
 *  the floats array to be 255.
 */
void write_pgm_normalize_float(const char *fname, const float *data,
                               int ncol, int nrow)
{
    FILE * f;
    int i,j;
    float v,max,min;

    /* check min and max values. If min is big than zero keep zero as min*/
    max = min = 0;
    for(i=0; i< nrow; i++)
        for(j=0; j<ncol; j++)
        {
            v = data[ j + i * ncol];
            if( v > max ) max = v;
            if( v < min ) min = v;
        }

    if( min < 0 ) fprintf(stderr,
        "Warning: negative values in '%s' are truncated to zero.\n",fname);

    /* open file */
    if( strcmp(fname,"-") == 0 ) f = stdout;
    else f = fopen(fname,"w");
    if( f == NULL ) errorf("Error: unable to open output image file ",fname);

    /* write header */
    fprintf(f,"P5\n");
    fprintf(f,"%u %u\n",ncol,nrow);
    fprintf(f,"%d\n",255);

    /* write data */
    for(i=0; i<nrow; i++)
        for(j=0; j<ncol; j++)
            fputc((unsigned char) (data[j+i*ncol]>0)
                  ? 255.0/max*data[j+i*ncol] : 0 ,f);

    /* close file if needed */
    if( f != stdout && fclose(f) == EOF )
        errorf("Error: unable to close file while writing PGM file ",fname);
}


/*---------------------------------------------------------------------------*/
/** Write an array of floats into a 3channel-8bit PPM image file. Normalize
 *   the maximum value of the floats array to be 255.
 */
void write_ppm_normalize_float(const char *fname, const float *rdata,
                               const float *gdata, const float *bdata,
                               int ncol, int nrow)
{
    FILE * f;
    int i,j;
    float v,max,min;

    /* check global min and max values in r,g,b channels. If min is big than
     *zero keep zero as min*/
    max = min = 0;

    /*r*/
    for(i=0; i< nrow; i++)
        for(j=0; j<ncol; j++)
        {
            v = rdata[ j + i * ncol];
            if( v > max ) max = v;
            if( v < min ) min = v;
        }

    /*g*/
    for(i=0; i< nrow; i++)
        for(j=0; j<ncol; j++)
        {
            v = gdata[ j + i * ncol];
            if( v > max ) max = v;
            if( v < min ) min = v;
        }

    /*b*/
    for(i=0; i< nrow; i++)
        for(j=0; j<ncol; j++)
        {
            v = bdata[ j + i * ncol];
            if( v > max ) max = v;
            if( v < min ) min = v;
        }

    if( min < 0 ) fprintf(stderr,
                 "Warning: negative values in '%s' are truncated to zero.\n",
                          fname);

    /* open file */
    if( strcmp(fname,"-") == 0 ) f = stdout;
    else f = fopen(fname,"w");
    if( f == NULL ) errorf("Error: unable to open output image file ",fname);

    /* write header */
    fprintf(f,"P6\n");
    fprintf(f,"%u %u\n",ncol,nrow);
    fprintf(f,"%d\n",255);

    /* write data */
    for(i=0; i<nrow; i++)
        for(j=0; j<ncol; j++)
        {
            static unsigned char color[3];
            /*red*/
            color[0] =(rdata[j+i*ncol]>0)? 255.0/max*rdata[j+i*ncol]:0;
            /*green*/
            color[1] =(gdata[j+i*ncol]>0)? 255.0/max*gdata[j+i*ncol]:0;
            /*blue*/
            color[2] =(bdata[j+i*ncol]>0)? 255.0/max*bdata[j+i*ncol]:0;
            (void) fwrite(color, 1, 3, f);
        }
    /* close file if needed */
    if( f != stdout && fclose(f) == EOF )
        errorf("Error: unable to close file while writing PPM file ",fname);
}
