#----------------------------------------------------------------------------
#
# "Point Spread Function Estimation from a Random Target"
#
# Copyright 2010-2011 mauricio delbracio (mdelbra@gmail.com)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# ---------------------------------------------------------------------------

# Makefile
# Non-parametric sub-pixel local point spread function estimation
# Author: Mauricio Delbracio <mdelbra@gmail.com>

# lapack and fftw libs/includes files should
# be in the following directories:
LIB_DIR = /opt/local/lib
LIB_INCLUDE_DIR = /opt/local/include

CSRC	= detect_pattern.c psf_estim.c psf_estim_main.c thin_plates.c lsd.c \
	image.c nnls.c io_pgm.c
OBJ	= $(CSRC:.c=.o)
BIN	= psf_estim

CFLAGS	= -ansi -pedantic -Wall -Wextra -Werror $(COPT)

default: $(BIN)

%.o	: %.c
	$(CC) $(CFLAGS) -I$(LIB_INCLUDE_DIR) -c -o $@ $<

$(BIN)  : io_pgm.o image.o lsd.o thin_plates.o detect_pattern.o nnls.o \
	 psf_estim.o psf_estim_main.o
	$(CC) $(CFLAGS) -L$(LIB_DIR) -o $@  $^ -lfftw3f -lm  -lblas -llapack

clean	:
	$(RM) $(OBJ)

doc	:
	doxygen doxygen.config

tar:
	tar -zcvf ../psfestim_1.3.tar.gz  --exclude=.* -C ../ psfestim_1.3/
