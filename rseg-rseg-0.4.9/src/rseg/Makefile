#    Makefile from common code directory of RSEG source tree
#
#    Copyright (C) 2011 University of Southern California and
#                       Andrew D. Smith
#
#    Authors: Qiang Song and Andrew D. Smith
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

SMITHLAB_CPP = ../smithlab_cpp
ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

COMMON_LIB = ../common
GSL_LIB = ../gsl
SAMTOOL_LIB = ../samtools

SOURCES = $(wildcard *.cpp)
OBJECTS = $(patsubst %.cpp,%.o,$(SOURCES))
PROGS =  rseg rseg-diff deadzones

INCLUDEDIRS = $(SMITHLAB_CPP) $(COMMON_LIB) $(GSL_LIB) $(SAMTOOL_LIB)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS = -lz -lm

CXX = g++
CC = gcc
CCFLAGS = -Wall -fPIC
CXXFLAGS = -Wall -fPIC -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g -lefence -L$(HOME)/lib

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

ifeq "$(shell uname)" "Darwin"
CXXFLAGS += -arch x86_64
CCFLAGS += -arch x86_64 -std=gnu89
ifeq "$(shell if [ `sysctl -n kern.osrelease | cut -d . -f 1` -ge 13 ];\
                  then echo 'true'; fi)" "true"
CXXFLAGS += -stdlib=libstdc++
endif
endif

all: $(PROGS)

rseg: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o) \
	$(addprefix $(COMMON_LIB)/, numerical_utils.o rseg_utils.o \
	Distro.o SplitDistro.o Smoothing.o BAMFile.o \
	SortGenomicRegion.o LoadReadsByRegion.o SelectBinSize.o ReadCounts.o \
	TwoStateScaleHMM.o TwoStateScaleResolveMixture.o \
	EvaluateBoundaries.o) \
	$(addprefix $(GSL_LIB)/, airy.o bessel_amp_phase.o bessel.o bessel_I0.o \
	bessel_I1.o bessel_In.o bessel_temme.o coerce.o elementary.o error.o exp.o \
	fdiv.o gamma.o gauss.o infnan.o log.o math.o mean.o message.o\
	mt.o poch.o pow_int.o psi.o rng.o stream.o trig.o variance.o \
	wmean.o wvariance.o zeta.o) \
	$(addprefix $(SAMTOOL_LIB)/, bgzf.o kstring.o bam_aux.o bam.o \
	bam_import.o sam.o bam_index.o bam_pileup.o bam_lpileup.o bam_md.o \
	razf.o faidx.o bedidx.o  bam_sort.o sam_header.o bam_reheader.o \
	kprobaln.o bam_cat.o  knetfile.o)


rseg-diff: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o RNG.o) \
	$(addprefix $(COMMON_LIB)/, numerical_utils.o rseg_utils.o \
	Distro.o SplitDistro.o Smoothing.o BAMFile.o \
	SortGenomicRegion.o LoadReadsByRegion.o SelectBinSize.o ReadCounts.o \
	TwoStateScaleSplitHMM.o TwoStateScaleSplitResolveMixture.o \
	ThreeStateScaleSplitHMM.o ThreeStateScaleSplitResolveMixture.o \
	EvaluateBoundaries.o) \
	$(addprefix $(GSL_LIB)/, airy.o bessel_amp_phase.o bessel.o bessel_I0.o \
	bessel_I1.o bessel_In.o bessel_temme.o coerce.o elementary.o error.o exp.o \
	fdiv.o gamma.o gauss.o infnan.o log.o math.o mean.o message.o\
	mt.o poch.o pow_int.o psi.o rng.o stream.o trig.o variance.o \
	wmean.o wvariance.o zeta.o) \
	$(addprefix $(SAMTOOL_LIB)/, bgzf.o kstring.o bam_aux.o bam.o \
	bam_import.o sam.o bam_index.o bam_pileup.o bam_lpileup.o bam_md.o \
	razf.o faidx.o bedidx.o  bam_sort.o sam_header.o bam_reheader.o \
	kprobaln.o bam_cat.o  knetfile.o)

deadzones: $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o smithlab_os.o \
	smithlab_utils.o OptionParser.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%.o: %.c
	$(CC) $(CCFLAGS) -c -o $@ $^ $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)


install: $(PROGS)
	@mkdir -p $(RSEG_ROOT)/bin
	@install -m 755 $(PROGS) $(RSEG_ROOT)/bin

clean:
	@-rm -f $(PROGS) *.o *~

.PHONY: clean
