
################################################################################
#
# Makefile for irs library blk
#

all:
	@ echo " "
	@ echo "Error:"
	@ echo "This Makefile must be called by the driving Makefile."
	@ echo " "

################################################################################

# Must use the bourne shell
#---------------------------
SHELL = /bin/sh

# Harmless suffix rule:
#----------------------
.SUFFIXES : .U

# library relative path name
#----------------------------
B = ../$(LIB_DIR)/libblk.a
BI = ../$(LIB_DIR)/libblk_inst.a

# Relative path to the directory of sources for the debugger
#-----------------------------------------------------------
DB  = ../../debug

# Relative path to the instrumented source directory from the source directory
#-----------------------------------------------------------------------------
PREDIR = ../../mccabe/pre
INSTDIR = ../../mccabe/inst

# Relative path to the instrumented source directory from the debug directory
#-----------------------------------------------------------------------------
INSTDIR_D = ../mccabe/inst

# Relative path to the source directory from the debug directory
#----------------------------------------------------------------
SRCDIR = ../sources/blk

# Relative path to the output file for lint
#------------------------------------------
# LO = ../../build/lint.out

# library objects list
#---------------------

libblk:$(B)(Blocks.o) $(B)(blkdiv.o) $(B)(blkdivin.o) $(B)(blkgen.o)              \
	  $(B)(blockin0.o) $(B)(domain.o) $(B)(domaincopy.o) $(B)(extents.o)      \
	  $(B)(initblk.o) $(B)(reblkdiv.o)
	  @if test -f object_file_flag; then  \
		(echo $$irs_ar ar $(B) *.o; \
		$$irs_ar $(B) *.o; \
		rm *.o object_file_flag); \
	  fi
	  @echo 'blk library current'

# creation deps
#---------------
$(B)(Blocks.o):Blocks.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/Blocks.c $(DB)/Blocks.c 

$(B)(blkdiv.o):blkdiv.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/blkdiv.c $(DB)/blkdiv.c 

$(B)(blkdivin.o):blkdivin.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/blkdivin.c $(DB)/blkdivin.c 

$(B)(blkgen.o):blkgen.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/blkgen.c $(DB)/blkgen.c 

$(B)(blockin0.o):blockin0.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/blockin0.c $(DB)/blockin0.c 

$(B)(domain.o):domain.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/domain.c $(DB)/domain.c 

$(B)(domaincopy.o):domaincopy.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/domaincopy.c $(DB)/domaincopy.c 

$(B)(extents.o):extents.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/extents.c $(DB)/extents.c 

$(B)(initblk.o):initblk.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/initblk.c $(DB)/initblk.c 

$(B)(reblkdiv.o):reblkdiv.c 
	$$irs_cc $<
	@touch object_file_flag
	@ln -sf $(SRCDIR)/reblkdiv.c $(DB)/reblkdiv.c 

#---------------------
# lint pass one target 
#---------------------
lint:Blocks.ln blkdiv.ln blkdivin.ln blkgen.ln blockin0.ln domain.ln              \
	  domaincopy.ln extents.ln initblk.ln reblkdiv.ln

#------------------------------
# individual lint .ln file deps
#------------------------------
Blocks.ln:Blocks.c ../../includes/irs.h ../../includes/irsblk.h                   \
	  ../../includes/irsbndry.h ../../includes/irsctl.h                       \
	  ../../includes/irsobjects.h ../../includes/irscom.h                     \
	  ../../includes/BlockMapping.h ../../includes/Extents.h                  \
	  ../../includes/FunctionTimer.h ../../includes/UserList.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON Blocks.c ***' >> $(LO)
	$$irs_lint Blocks.c 2>&1 | tee -a $(LO)

blkdiv.ln:blkdiv.c ../../includes/irs.h ../../includes/irsgen.h                   \
	  ../../includes/irsblk.h ../../includes/irsreg.h                         \
	  ../../includes/BlockMapping.h ../../includes/Extents.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON blkdiv.c ***' >> $(LO)
	$$irs_lint blkdiv.c 2>&1 | tee -a $(LO)

blkdivin.ln:blkdivin.c ../../includes/irs.h ../../includes/irsblk.h               \
	  ../../includes/irscom.h ../../includes/irsctl.h                         \
	  ../../includes/BlockMapping.h ../../includes/Extents.h                  \
	  ../../includes/FunctionTimer.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON blkdivin.c ***' >> $(LO)
	$$irs_lint blkdivin.c 2>&1 | tee -a $(LO)

blkgen.ln:blkgen.c ../../includes/irs.h ../../includes/irssrc.h                   \
	  ../../includes/irsbndry.h ../../includes/irsblk.h                       \
	  ../../includes/irsgen.h ../../includes/irsreg.h                         \
	  ../../includes/irsrgst.h ../../includes/BlockMapping.h                  \
	  ../../includes/Command.h ../../includes/Source.h                        \
	  ../../includes/Domain.h ../../includes/Rgst.h ../../includes/Region.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON blkgen.c ***' >> $(LO)
	$$irs_lint blkgen.c 2>&1 | tee -a $(LO)

blockin0.ln:blockin0.c ../../includes/irs.h ../../includes/irsblk.h               \
	  ../../includes/irsbndry.h ../../includes/irscom.h                       \
	  ../../includes/irsctl.h ../../includes/BlockMapping.h                   \
	  ../../includes/FunctionTimer.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON blockin0.c ***' >> $(LO)
	$$irs_lint blockin0.c 2>&1 | tee -a $(LO)

domain.ln:domain.c ../../includes/irs.h ../../includes/irsblk.h                   \
	  ../../includes/irscom.h ../../includes/irsobjects.h                     \
	  ../../includes/irsparm.h ../../includes/irsreg.h                        \
	  ../../includes/irsrgst.h ../../includes/BoundaryCondition.h             \
	  ../../includes/Boundary.h ../../includes/Command.h                      \
	  ../../includes/Domain.h ../../includes/Hash.h                           \
	  ../../includes/ProblemArray.h ../../includes/FunctionTimer.h            \
	  ../../includes/Rgst.h ../../includes/SpeciesFraction.h                  \
	  ../../includes/TimeStepControl.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON domain.c ***' >> $(LO)
	$$irs_lint domain.c 2>&1 | tee -a $(LO)

domaincopy.ln:domaincopy.c ../../includes/irs.h ../../includes/irssrc.h           \
	  ../../includes/irsblk.h ../../includes/irsbndry.h                       \
	  ../../includes/irsrgst.h ../../includes/irsreg.h ../../includes/pint.h  \
	  ../../includes/ProblemArray.h ../../includes/Rgst.h                     \
	  ../../includes/Region.h ../../includes/SpeciesFraction.h                \
	  ../../includes/Extents.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON domaincopy.c ***' >> $(LO)
	$$irs_lint domaincopy.c 2>&1 | tee -a $(LO)

extents.ln:extents.c ../../includes/irs.h ../../includes/irsblk.h                 \
	  ../../includes/Extents.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON extents.c ***' >> $(LO)
	$$irs_lint extents.c 2>&1 | tee -a $(LO)

initblk.ln:initblk.c ../../includes/irs.h ../../includes/irsgen.h                 \
	  ../../includes/irsblk.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON initblk.c ***' >> $(LO)
	$$irs_lint initblk.c 2>&1 | tee -a $(LO)

reblkdiv.ln:reblkdiv.c ../../includes/irs.h ../../includes/irsblk.h               \
	  ../../includes/irscom.h ../../includes/irsdmp.h ../../includes/irsgen.h \
	  ../../includes/irsrgst.h ../../includes/BlockMapping.h                  \
	  ../../includes/Extents.h ../../includes/NodeList.h
	  @echo '  *** RUNNING $(SYNTAX_CHECKER) ON reblkdiv.c ***' >> $(LO)
	$$irs_lint reblkdiv.c 2>&1 | tee -a $(LO)

#------------------------------------------
# preprocess files and create .i files 
#------------------------------------------
preprocess:$(PREDIR)/Blocks.inst.i $(PREDIR)/blkdiv.inst.i                        \
	  $(PREDIR)/blkdivin.inst.i $(PREDIR)/blkgen.inst.i                       \
	  $(PREDIR)/blockin0.inst.i $(PREDIR)/domain.inst.i                       \
	  $(PREDIR)/domaincopy.inst.i $(PREDIR)/extents.inst.i                    \
	  $(PREDIR)/initblk.inst.i $(PREDIR)/reblkdiv.inst.i

#---------------------------------------
# individual preprocessed .i file deps
#---------------------------------------
$(PREDIR)/Blocks.inst.i:Blocks.c ../../includes/irs.h ../../includes/irsblk.h     \
	  ../../includes/irsbndry.h ../../includes/irsctl.h                       \
	  ../../includes/irsobjects.h ../../includes/irscom.h                     \
	  ../../includes/BlockMapping.h ../../includes/Extents.h                  \
	  ../../includes/FunctionTimer.h ../../includes/UserList.h
	@rm -f Blocks.i
	$$irs_pre Blocks.c >> Blocks.i
	@mv -f Blocks.i $(PREDIR)/Blocks.inst.i

$(PREDIR)/blkdiv.inst.i:blkdiv.c ../../includes/irs.h ../../includes/irsgen.h     \
	  ../../includes/irsblk.h ../../includes/irsreg.h                         \
	  ../../includes/BlockMapping.h ../../includes/Extents.h
	@rm -f blkdiv.i
	$$irs_pre blkdiv.c >> blkdiv.i
	@mv -f blkdiv.i $(PREDIR)/blkdiv.inst.i

$(PREDIR)/blkdivin.inst.i:blkdivin.c ../../includes/irs.h ../../includes/irsblk.h \
	  ../../includes/irscom.h ../../includes/irsctl.h                         \
	  ../../includes/BlockMapping.h ../../includes/Extents.h                  \
	  ../../includes/FunctionTimer.h
	@rm -f blkdivin.i
	$$irs_pre blkdivin.c >> blkdivin.i
	@mv -f blkdivin.i $(PREDIR)/blkdivin.inst.i

$(PREDIR)/blkgen.inst.i:blkgen.c ../../includes/irs.h ../../includes/irssrc.h     \
	  ../../includes/irsbndry.h ../../includes/irsblk.h                       \
	  ../../includes/irsgen.h ../../includes/irsreg.h                         \
	  ../../includes/irsrgst.h ../../includes/BlockMapping.h                  \
	  ../../includes/Command.h ../../includes/Source.h                        \
	  ../../includes/Domain.h ../../includes/Rgst.h ../../includes/Region.h
	@rm -f blkgen.i
	$$irs_pre blkgen.c >> blkgen.i
	@mv -f blkgen.i $(PREDIR)/blkgen.inst.i

$(PREDIR)/blockin0.inst.i:blockin0.c ../../includes/irs.h ../../includes/irsblk.h \
	  ../../includes/irsbndry.h ../../includes/irscom.h                       \
	  ../../includes/irsctl.h ../../includes/BlockMapping.h                   \
	  ../../includes/FunctionTimer.h
	@rm -f blockin0.i
	$$irs_pre blockin0.c >> blockin0.i
	@mv -f blockin0.i $(PREDIR)/blockin0.inst.i

$(PREDIR)/domain.inst.i:domain.c ../../includes/irs.h ../../includes/irsblk.h     \
	  ../../includes/irscom.h ../../includes/irsobjects.h                     \
	  ../../includes/irsparm.h ../../includes/irsreg.h                        \
	  ../../includes/irsrgst.h ../../includes/BoundaryCondition.h             \
	  ../../includes/Boundary.h ../../includes/Command.h                      \
	  ../../includes/Domain.h ../../includes/Hash.h                           \
	  ../../includes/ProblemArray.h ../../includes/FunctionTimer.h            \
	  ../../includes/Rgst.h ../../includes/SpeciesFraction.h                  \
	  ../../includes/TimeStepControl.h
	@rm -f domain.i
	$$irs_pre domain.c >> domain.i
	@mv -f domain.i $(PREDIR)/domain.inst.i

$(PREDIR)/domaincopy.inst.i:domaincopy.c ../../includes/irs.h                     \
	  ../../includes/irssrc.h ../../includes/irsblk.h                         \
	  ../../includes/irsbndry.h ../../includes/irsrgst.h                      \
	  ../../includes/irsreg.h ../../includes/pint.h                           \
	  ../../includes/ProblemArray.h ../../includes/Rgst.h                     \
	  ../../includes/Region.h ../../includes/SpeciesFraction.h                \
	  ../../includes/Extents.h
	@rm -f domaincopy.i
	$$irs_pre domaincopy.c >> domaincopy.i
	@mv -f domaincopy.i $(PREDIR)/domaincopy.inst.i

$(PREDIR)/extents.inst.i:extents.c ../../includes/irs.h ../../includes/irsblk.h   \
	  ../../includes/Extents.h
	@rm -f extents.i
	$$irs_pre extents.c >> extents.i
	@mv -f extents.i $(PREDIR)/extents.inst.i

$(PREDIR)/initblk.inst.i:initblk.c ../../includes/irs.h ../../includes/irsgen.h   \
	  ../../includes/irsblk.h
	@rm -f initblk.i
	$$irs_pre initblk.c >> initblk.i
	@mv -f initblk.i $(PREDIR)/initblk.inst.i

$(PREDIR)/reblkdiv.inst.i:reblkdiv.c ../../includes/irs.h ../../includes/irsblk.h \
	  ../../includes/irscom.h ../../includes/irsdmp.h ../../includes/irsgen.h \
	  ../../includes/irsrgst.h ../../includes/BlockMapping.h                  \
	  ../../includes/Extents.h ../../includes/NodeList.h
	@rm -f reblkdiv.i
	$$irs_pre reblkdiv.c >> reblkdiv.i
	@mv -f reblkdiv.i $(PREDIR)/reblkdiv.inst.i

#-----------------------------------------------------------
# post-process files and create library from .i files 
#-----------------------------------------------------------
postprocess:libblk_inst

libblk_inst:$(BI)(Blocks.inst.o) $(BI)(blkdiv.inst.o) $(BI)(blkdivin.inst.o)      \
	  $(BI)(blkgen.inst.o) $(BI)(blockin0.inst.o) $(BI)(domain.inst.o)        \
	  $(BI)(domaincopy.inst.o) $(BI)(extents.inst.o) $(BI)(initblk.inst.o)    \
	  $(BI)(reblkdiv.inst.o)
	  @if test -f object_file_flag; then  \
		(echo $$irs_ar ar $(BI) *.o; \
		$$irs_ar $(BI) *.o; \
		rm *.o object_file_flag; \
		echo ln -sf libblk_inst.a $(B); \
		ln -sf libblk_inst.a $(B)); \
	  fi
	  @echo 'libblk_inst library current'


#-------------------------------------------------------------------------------
# individual rules to compile .i files into .o files within instrumented library
#-------------------------------------------------------------------------------
$(BI)(Blocks.inst.o):$(INSTDIR)/Blocks.inst.i
	$$irs_post $(INSTDIR)/Blocks.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/Blocks.inst.i $(DB)/Blocks.c 

$(BI)(blkdiv.inst.o):$(INSTDIR)/blkdiv.inst.i
	$$irs_post $(INSTDIR)/blkdiv.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/blkdiv.inst.i $(DB)/blkdiv.c 

$(BI)(blkdivin.inst.o):$(INSTDIR)/blkdivin.inst.i
	$$irs_post $(INSTDIR)/blkdivin.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/blkdivin.inst.i $(DB)/blkdivin.c 

$(BI)(blkgen.inst.o):$(INSTDIR)/blkgen.inst.i
	$$irs_post $(INSTDIR)/blkgen.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/blkgen.inst.i $(DB)/blkgen.c 

$(BI)(blockin0.inst.o):$(INSTDIR)/blockin0.inst.i
	$$irs_post $(INSTDIR)/blockin0.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/blockin0.inst.i $(DB)/blockin0.c 

$(BI)(domain.inst.o):$(INSTDIR)/domain.inst.i
	$$irs_post $(INSTDIR)/domain.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/domain.inst.i $(DB)/domain.c 

$(BI)(domaincopy.inst.o):$(INSTDIR)/domaincopy.inst.i
	$$irs_post $(INSTDIR)/domaincopy.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/domaincopy.inst.i $(DB)/domaincopy.c 

$(BI)(extents.inst.o):$(INSTDIR)/extents.inst.i
	$$irs_post $(INSTDIR)/extents.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/extents.inst.i $(DB)/extents.c 

$(BI)(initblk.inst.o):$(INSTDIR)/initblk.inst.i
	$$irs_post $(INSTDIR)/initblk.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/initblk.inst.i $(DB)/initblk.c 

$(BI)(reblkdiv.inst.o):$(INSTDIR)/reblkdiv.inst.i
	$$irs_post $(INSTDIR)/reblkdiv.inst.i
	@touch object_file_flag
	@ln -sf $(INSTDIR_D)/reblkdiv.inst.i $(DB)/reblkdiv.c 

DEP_SRCS = \
	Blocks.c blkdiv.c blkdivin.c blkgen.c blockin0.c domain.c domaincopy.c  \
	extents.c initblk.c reblkdiv.c  

depend:
	@cp -f Makefile.blk Makefile
	makedepend -- $(DEPFLAGS) -- $(DEP_SRCS)
	irs_fixdepend Makefile

# DO NOT DELETE

Blocks.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
Blocks.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
Blocks.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
Blocks.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
Blocks.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
Blocks.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
Blocks.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
Blocks.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
Blocks.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
Blocks.o: /usr/include/sys/types.h /usr/include/time.h
Blocks.o: /usr/include/sys/select.h /usr/include/bits/select.h
Blocks.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
Blocks.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
Blocks.o: /usr/include/alloca.h /usr/include/sys/times.h
Blocks.o: ../../includes/irsdefs.h ../../includes/irsctl.h
Blocks.o: ../../includes/Rgst.h ../../includes/Hash.h ../../includes/Domain.h
Blocks.o: ../../includes/pint.h ../../includes/SpeciesFraction.h
Blocks.o: ../../includes/Source.h ../../includes/BoundaryCondition.h
Blocks.o: ../../includes/Boundary.h ../../includes/mach.h
Blocks.o: /usr/include/ctype.h /usr/include/xlocale.h /usr/include/math.h
Blocks.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
Blocks.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
Blocks.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
Blocks.o: /usr/include/bits/mathcalls.h /usr/include/string.h
Blocks.o: /usr/include/malloc.h ../../includes/irsblk.h
Blocks.o: ../../includes/Extents.h ../../includes/irsbndry.h
Blocks.o: ../../includes/irsobjects.h ../../includes/irscom.h
Blocks.o: ../../includes/Combuffer.h ../../includes/FunctionTimer.h
Blocks.o: ../../includes/NodeList.h ../../includes/BlockMapping.h
Blocks.o: ../../includes/BndMapping.h ../../includes/UserList.h
blkdiv.o: /usr/include/unistd.h /usr/include/features.h
blkdiv.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
blkdiv.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
blkdiv.o: /usr/include/bits/posix_opt.h /usr/include/bits/environments.h
blkdiv.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
blkdiv.o: /usr/include/bits/confname.h /usr/include/getopt.h
blkdiv.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/libio.h
blkdiv.o: /usr/include/_G_config.h /usr/include/wchar.h
blkdiv.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blkdiv.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
blkdiv.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
blkdiv.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
blkdiv.o: /usr/include/sys/types.h /usr/include/time.h
blkdiv.o: /usr/include/sys/select.h /usr/include/bits/select.h
blkdiv.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
blkdiv.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
blkdiv.o: /usr/include/alloca.h /usr/include/sys/times.h
blkdiv.o: ../../includes/irsdefs.h ../../includes/irsctl.h
blkdiv.o: ../../includes/Rgst.h ../../includes/Hash.h ../../includes/Domain.h
blkdiv.o: ../../includes/pint.h ../../includes/SpeciesFraction.h
blkdiv.o: ../../includes/Source.h ../../includes/BoundaryCondition.h
blkdiv.o: ../../includes/Boundary.h ../../includes/mach.h
blkdiv.o: /usr/include/ctype.h /usr/include/xlocale.h /usr/include/math.h
blkdiv.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
blkdiv.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
blkdiv.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
blkdiv.o: /usr/include/bits/mathcalls.h /usr/include/string.h
blkdiv.o: /usr/include/malloc.h ../../includes/irsgen.h
blkdiv.o: ../../includes/NodeList.h ../../includes/Extents.h
blkdiv.o: ../../includes/irsblk.h ../../includes/irsreg.h
blkdiv.o: ../../includes/Region.h ../../includes/BlockMapping.h
blkdiv.o: ../../includes/BndMapping.h
blkdivin.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
blkdivin.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
blkdivin.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
blkdivin.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
blkdivin.o: /usr/include/libio.h /usr/include/_G_config.h
blkdivin.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
blkdivin.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
blkdivin.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
blkdivin.o: /usr/include/endian.h /usr/include/bits/endian.h
blkdivin.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
blkdivin.o: /usr/include/time.h /usr/include/sys/select.h
blkdivin.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
blkdivin.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
blkdivin.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
blkdivin.o: /usr/include/sys/times.h ../../includes/irsdefs.h
blkdivin.o: ../../includes/irsctl.h ../../includes/Rgst.h
blkdivin.o: ../../includes/Hash.h ../../includes/Domain.h
blkdivin.o: ../../includes/pint.h ../../includes/SpeciesFraction.h
blkdivin.o: ../../includes/Source.h ../../includes/BoundaryCondition.h
blkdivin.o: ../../includes/Boundary.h ../../includes/mach.h
blkdivin.o: /usr/include/ctype.h /usr/include/xlocale.h /usr/include/math.h
blkdivin.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
blkdivin.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
blkdivin.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
blkdivin.o: /usr/include/bits/mathcalls.h /usr/include/string.h
blkdivin.o: /usr/include/malloc.h ../../includes/irsblk.h
blkdivin.o: ../../includes/Extents.h ../../includes/irscom.h
blkdivin.o: ../../includes/Combuffer.h ../../includes/FunctionTimer.h
blkdivin.o: ../../includes/NodeList.h ../../includes/BlockMapping.h
blkdivin.o: ../../includes/BndMapping.h
blkgen.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
blkgen.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
blkgen.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
blkgen.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
blkgen.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
blkgen.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
blkgen.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
blkgen.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
blkgen.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
blkgen.o: /usr/include/sys/types.h /usr/include/time.h
blkgen.o: /usr/include/sys/select.h /usr/include/bits/select.h
blkgen.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
blkgen.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
blkgen.o: /usr/include/alloca.h /usr/include/sys/times.h
blkgen.o: ../../includes/irsdefs.h ../../includes/irsctl.h
blkgen.o: ../../includes/Rgst.h ../../includes/Hash.h ../../includes/Domain.h
blkgen.o: ../../includes/pint.h ../../includes/SpeciesFraction.h
blkgen.o: ../../includes/Source.h ../../includes/BoundaryCondition.h
blkgen.o: ../../includes/Boundary.h ../../includes/mach.h
blkgen.o: /usr/include/ctype.h /usr/include/xlocale.h /usr/include/math.h
blkgen.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
blkgen.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
blkgen.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
blkgen.o: /usr/include/bits/mathcalls.h /usr/include/string.h
blkgen.o: /usr/include/malloc.h ../../includes/irssrc.h
blkgen.o: ../../includes/NodeList.h ../../includes/Extents.h
blkgen.o: ../../includes/irsbndry.h ../../includes/irsblk.h
blkgen.o: ../../includes/irsgen.h ../../includes/irsreg.h
blkgen.o: ../../includes/Region.h ../../includes/irsrgst.h
blkgen.o: ../../includes/BlockMapping.h ../../includes/BndMapping.h
blkgen.o: ../../includes/Command.h
blockin0.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
blockin0.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
blockin0.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
blockin0.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
blockin0.o: /usr/include/libio.h /usr/include/_G_config.h
blockin0.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
blockin0.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
blockin0.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
blockin0.o: /usr/include/endian.h /usr/include/bits/endian.h
blockin0.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
blockin0.o: /usr/include/time.h /usr/include/sys/select.h
blockin0.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
blockin0.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
blockin0.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
blockin0.o: /usr/include/sys/times.h ../../includes/irsdefs.h
blockin0.o: ../../includes/irsctl.h ../../includes/Rgst.h
blockin0.o: ../../includes/Hash.h ../../includes/Domain.h
blockin0.o: ../../includes/pint.h ../../includes/SpeciesFraction.h
blockin0.o: ../../includes/Source.h ../../includes/BoundaryCondition.h
blockin0.o: ../../includes/Boundary.h ../../includes/mach.h
blockin0.o: /usr/include/ctype.h /usr/include/xlocale.h /usr/include/math.h
blockin0.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
blockin0.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
blockin0.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
blockin0.o: /usr/include/bits/mathcalls.h /usr/include/string.h
blockin0.o: /usr/include/malloc.h ../../includes/irsblk.h
blockin0.o: ../../includes/Extents.h ../../includes/irsbndry.h
blockin0.o: ../../includes/irscom.h ../../includes/Combuffer.h
blockin0.o: ../../includes/FunctionTimer.h ../../includes/NodeList.h
blockin0.o: ../../includes/BlockMapping.h ../../includes/BndMapping.h
domain.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
domain.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
domain.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
domain.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
domain.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
domain.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
domain.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
domain.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
domain.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
domain.o: /usr/include/sys/types.h /usr/include/time.h
domain.o: /usr/include/sys/select.h /usr/include/bits/select.h
domain.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
domain.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
domain.o: /usr/include/alloca.h /usr/include/sys/times.h
domain.o: ../../includes/irsdefs.h ../../includes/irsctl.h
domain.o: ../../includes/Rgst.h ../../includes/Hash.h ../../includes/Domain.h
domain.o: ../../includes/pint.h ../../includes/SpeciesFraction.h
domain.o: ../../includes/Source.h ../../includes/BoundaryCondition.h
domain.o: ../../includes/Boundary.h ../../includes/mach.h
domain.o: /usr/include/ctype.h /usr/include/xlocale.h /usr/include/math.h
domain.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
domain.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
domain.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
domain.o: /usr/include/bits/mathcalls.h /usr/include/string.h
domain.o: /usr/include/malloc.h ../../includes/irsblk.h
domain.o: ../../includes/Extents.h ../../includes/irscom.h
domain.o: ../../includes/Combuffer.h ../../includes/FunctionTimer.h
domain.o: ../../includes/NodeList.h ../../includes/irsobjects.h
domain.o: ../../includes/irsparm.h ../../includes/irsreg.h
domain.o: ../../includes/Region.h ../../includes/irsrgst.h
domain.o: ../../includes/Command.h ../../includes/ProblemArray.h
domain.o: ../../includes/TimeStepControl.h
domaincopy.o: ../../includes/irs.h /usr/include/stdio.h
domaincopy.o: /usr/include/features.h /usr/include/sys/cdefs.h
domaincopy.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
domaincopy.o: /usr/include/gnu/stubs-64.h /usr/include/bits/types.h
domaincopy.o: /usr/include/bits/typesizes.h /usr/include/libio.h
domaincopy.o: /usr/include/_G_config.h /usr/include/wchar.h
domaincopy.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
domaincopy.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
domaincopy.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
domaincopy.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
domaincopy.o: /usr/include/sys/types.h /usr/include/time.h
domaincopy.o: /usr/include/sys/select.h /usr/include/bits/select.h
domaincopy.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
domaincopy.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
domaincopy.o: /usr/include/alloca.h /usr/include/sys/times.h
domaincopy.o: ../../includes/irsdefs.h ../../includes/irsctl.h
domaincopy.o: ../../includes/Rgst.h ../../includes/Hash.h
domaincopy.o: ../../includes/Domain.h ../../includes/pint.h
domaincopy.o: ../../includes/SpeciesFraction.h ../../includes/Source.h
domaincopy.o: ../../includes/BoundaryCondition.h ../../includes/Boundary.h
domaincopy.o: ../../includes/mach.h /usr/include/ctype.h
domaincopy.o: /usr/include/xlocale.h /usr/include/math.h
domaincopy.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
domaincopy.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
domaincopy.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
domaincopy.o: /usr/include/bits/mathcalls.h /usr/include/string.h
domaincopy.o: /usr/include/malloc.h ../../includes/irssrc.h
domaincopy.o: ../../includes/NodeList.h ../../includes/Extents.h
domaincopy.o: ../../includes/irsblk.h ../../includes/irsbndry.h
domaincopy.o: ../../includes/irsrgst.h ../../includes/irsreg.h
domaincopy.o: ../../includes/Region.h ../../includes/ProblemArray.h
extents.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
extents.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
extents.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
extents.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
extents.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
extents.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
extents.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
extents.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
extents.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
extents.o: /usr/include/sys/types.h /usr/include/time.h
extents.o: /usr/include/sys/select.h /usr/include/bits/select.h
extents.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
extents.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
extents.o: /usr/include/alloca.h /usr/include/sys/times.h
extents.o: ../../includes/irsdefs.h ../../includes/irsctl.h
extents.o: ../../includes/Rgst.h ../../includes/Hash.h
extents.o: ../../includes/Domain.h ../../includes/pint.h
extents.o: ../../includes/SpeciesFraction.h ../../includes/Source.h
extents.o: ../../includes/BoundaryCondition.h ../../includes/Boundary.h
extents.o: ../../includes/mach.h /usr/include/ctype.h /usr/include/xlocale.h
extents.o: /usr/include/math.h /usr/include/bits/huge_val.h
extents.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
extents.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
extents.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
extents.o: /usr/include/string.h /usr/include/malloc.h
extents.o: ../../includes/irsblk.h ../../includes/Extents.h
initblk.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
initblk.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
initblk.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
initblk.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
initblk.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
initblk.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
initblk.o: /usr/include/stdlib.h /usr/include/bits/waitflags.h
initblk.o: /usr/include/bits/waitstatus.h /usr/include/endian.h
initblk.o: /usr/include/bits/endian.h /usr/include/bits/byteswap.h
initblk.o: /usr/include/sys/types.h /usr/include/time.h
initblk.o: /usr/include/sys/select.h /usr/include/bits/select.h
initblk.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
initblk.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
initblk.o: /usr/include/alloca.h /usr/include/sys/times.h
initblk.o: ../../includes/irsdefs.h ../../includes/irsctl.h
initblk.o: ../../includes/Rgst.h ../../includes/Hash.h
initblk.o: ../../includes/Domain.h ../../includes/pint.h
initblk.o: ../../includes/SpeciesFraction.h ../../includes/Source.h
initblk.o: ../../includes/BoundaryCondition.h ../../includes/Boundary.h
initblk.o: ../../includes/mach.h /usr/include/ctype.h /usr/include/xlocale.h
initblk.o: /usr/include/math.h /usr/include/bits/huge_val.h
initblk.o: /usr/include/bits/huge_valf.h /usr/include/bits/huge_vall.h
initblk.o: /usr/include/bits/inf.h /usr/include/bits/nan.h
initblk.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
initblk.o: /usr/include/string.h /usr/include/malloc.h
initblk.o: ../../includes/irsgen.h ../../includes/NodeList.h
initblk.o: ../../includes/Extents.h ../../includes/irsblk.h
reblkdiv.o: ../../includes/irs.h /usr/include/stdio.h /usr/include/features.h
reblkdiv.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
reblkdiv.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-64.h
reblkdiv.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
reblkdiv.o: /usr/include/libio.h /usr/include/_G_config.h
reblkdiv.o: /usr/include/wchar.h /usr/include/bits/stdio_lim.h
reblkdiv.o: /usr/include/bits/sys_errlist.h /usr/include/stdlib.h
reblkdiv.o: /usr/include/bits/waitflags.h /usr/include/bits/waitstatus.h
reblkdiv.o: /usr/include/endian.h /usr/include/bits/endian.h
reblkdiv.o: /usr/include/bits/byteswap.h /usr/include/sys/types.h
reblkdiv.o: /usr/include/time.h /usr/include/sys/select.h
reblkdiv.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
reblkdiv.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
reblkdiv.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
reblkdiv.o: /usr/include/sys/times.h ../../includes/irsdefs.h
reblkdiv.o: ../../includes/irsctl.h ../../includes/Rgst.h
reblkdiv.o: ../../includes/Hash.h ../../includes/Domain.h
reblkdiv.o: ../../includes/pint.h ../../includes/SpeciesFraction.h
reblkdiv.o: ../../includes/Source.h ../../includes/BoundaryCondition.h
reblkdiv.o: ../../includes/Boundary.h ../../includes/mach.h
reblkdiv.o: /usr/include/ctype.h /usr/include/xlocale.h /usr/include/math.h
reblkdiv.o: /usr/include/bits/huge_val.h /usr/include/bits/huge_valf.h
reblkdiv.o: /usr/include/bits/huge_vall.h /usr/include/bits/inf.h
reblkdiv.o: /usr/include/bits/nan.h /usr/include/bits/mathdef.h
reblkdiv.o: /usr/include/bits/mathcalls.h /usr/include/string.h
reblkdiv.o: /usr/include/malloc.h ../../includes/irsblk.h
reblkdiv.o: ../../includes/Extents.h ../../includes/irscom.h
reblkdiv.o: ../../includes/Combuffer.h ../../includes/FunctionTimer.h
reblkdiv.o: ../../includes/NodeList.h ../../includes/irsdmp.h
reblkdiv.o: ../../includes/irsgen.h ../../includes/irsrgst.h
reblkdiv.o: ../../includes/BlockMapping.h ../../includes/BndMapping.h
