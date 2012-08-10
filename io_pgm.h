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


/**
 * @file io_pgm.h
 * @brief library header for basic image processing.
 * @author Mauricio Delbracio  (mdelbra@gmail.com)
 */

#ifndef IO_PGM_H_
#define IO_PGM_H_


float *read_pgm_float(const char * fname, int * ncol, int *nrow);

void write_pgm_float(const char *fname, const float *data, int ncol, int nrow);

void write_pgm_normalize_float(const char *fname, const float *data,
                               int ncol, int nrow);

void write_ppm_normalize_float(const char *fname, const float *rdata,
                               const float *gdata, const float *bdata,
                               int ncol, int nrow);


#endif
