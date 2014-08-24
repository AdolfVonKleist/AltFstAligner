# Makefile
#
# Copyright (c) [2012-], Josef Robert Novak
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
#  modification, are permitted #provided that the following conditions
#  are met:
#
#  * Redistributions of source code must retain the above copyright 
#    notice, this list of conditions and the following disclaimer.
#    * Redistributions in binary form must reproduce the above 
#    copyright notice, this list of #conditions and the following 
#    disclaimer in the documentation and/or other materials provided 
#    with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
# OF THE POSSIBILITY OF SUCH DAMAGE.
CC=g++

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
LIBS=-lfst /usr/local/lib/fst/libfstfar.dylib
endif

ifeq ($(UNAME), Linux)
LIBS=-lfst -ldl -lm -lrt -lfstfar
endif
CC=g++ -std=c++0x -Wall -Wno-sign-compare

OUT=$(OBJS)

CFLAGS ?= -O2     
THIRD_PARTIES_INCLUDE ?= -I3rdparty/utfcpp

EXTRA=$(CFLAGS) $(THIRD_PARTIES_INCLUDE) $(CXXFLAGS) $(CPPFLAGS) #-fopenmp

%.o: %.cc
	$(CC) $(EXTRA) -c $(<) -o $(@) 

all: altfst-align

altfst-align: Aligner.h DictionaryReader.h altfst-align.o
	$(CC) $(EXTRA) $(LIBS) -L/usr/local/lib/fst \
	-funroll-loops -I. -ffast-math altfst-align.o \
	-o altfst-align $(LIBS)

install: all
	install -c altfst-align ../altfst-align
clean:
	$(RM) $(OUT) $(TMP) *.o altfst-align ../altfst-align
