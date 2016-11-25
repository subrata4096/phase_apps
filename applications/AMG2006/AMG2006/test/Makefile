#BHEADER**********************************************************************
# Copyright (c) 2006   The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# Written by the HYPRE team. UCRL-CODE-222953.
# All rights reserved.
#
# This file is part of HYPRE (see http://www.llnl.gov/CASC/hypre/).
# Please see the COPYRIGHT_and_LICENSE file for the copyright notice, 
# disclaimer, contact information and the GNU Lesser General Public License.
#
# HYPRE is free software; you can redistribute it and/or modify it under the
# terms of the GNU General Public License (as published by the Free Software 
# Foundation) version 2.1 dated February 1999.
#
# HYPRE is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the terms and conditions of the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
# $Revision: 1.50 $
#EHEADER**********************************************************************

.SUFFIXES:
.SUFFIXES: .c .o
                                                                                               
srcdir = .
                                                                                               
include ../Makefile.include
                                                                                               
CFLAGS = \
 -I..\
 -I../utilities\
 -I../struct_mv\
 -I../sstruct_mv\
 -I../IJ_mv\
 -I../seq_mv\
 -I../parcsr_mv\
 -I../parcsr_ls\
 -I../krylov\
 ${INCLUDE_CFLAGS}\
 -DHYPRE_TIMING
                                                                                               
 #/g/g90/mitra3/work/for_cab/automaded_LA_profiling/install/lib/libstracker.a\
 #-L/g/g90/mitra3/work/for_cab/automaded_LA_24Oct/install/lib\
 #-lstracker
 #-L/collab/usr/global/tools/totalview/r/toolworks/totalview.8.11.0-3/linux-x86/lib\
 -ltvheap
 #-L/usr/global/tools/adept/callpath/chaos_5_x86_64_ib/lib
 #-lcallpath
 #-L/usr/global/tools/adept-utils/chaos_5_x86_64_ib/lib
 #-ladept_utils
LFLAGS =\
 -L.\
 -L../parcsr_ls\
 -L../parcsr_mv\
 -L../IJ_mv\
 -L../seq_mv\
 -L../sstruct_mv\
 -L../struct_mv\
 -L../krylov\
 -L../utilities\
 -lparcsr_ls\
 -lparcsr_mv\
 -lseq_mv\
 -lsstruct_mv\
 -lIJ_mv\
 -lHYPRE_struct_mv\
 -lkrylov\
 -lHYPRE_utilities\
 -rdynamic\
 #/g/g90/mitra3/work/code/automaded_profiler/install/lib/libstracker.a -L/usr/global/tools/adept/callpath/chaos_5_x86_64_ib/lib -lcallpath -L/usr/global/tools/adept-utils/chaos_5_x86_64_ib/lib -ladept_utils\
 ${INCLUDE_LFLAGS}

##################################################################
# Targets
##################################################################
                                                                                                                                         
all: amg2006
                                                                                                                                         
clean:
	@rm -f *.o
                                                                                                                                         
veryclean: clean
	@rm -f amg2006
                                                                                                                                         
##################################################################
# Rules
##################################################################
                                                                                                                                         
amg2006: amg2006.o
	@echo  "Linking" $@ "... "
	${CC} -g -O3 -o amg2006 amg2006.o ${LFLAGS}
                                                                                                                                         
##################################################################
# Generic rules
##################################################################
                                                                                                                                         
.c.o:
	${CC} -g -O3 -o $@ -c ${CFLAGS} $<
                                                               
