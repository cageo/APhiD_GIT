#######################################
#                                     #
#  Makefile for APhiD                 #
#                                     #
#  Chester J Weiss                    #
#  Virginia Tech Dept of Geosciences  #
#  cjweiss@vt.edu                     #
#  chester.weiss@gmail.com            #
#                                     #         
#######################################

###############################################################################
#  BUILD OPTIONS                                                              #
#                                                                             #
#  make serial (default) : install serial version in $(DEST)                  #
#  make parallel         : install serial version in $(DEST)                  #
#  make man              : install man page in $(MANPATH), 'man aphid' to see #
#  make clean            : delete *.o, *.mod files                            #
#                                                                             #
###############################################################################

###############################################################################
#  Set compute platform ( mac | altix | cabernet | hess | anodyne | ...  ? )  #
###############################################################################

PLATFORM = mac_gfortran

#######################################################################
#  Compiler, options, and directory info for platforms defined above  #
#######################################################################

ifeq ($(PLATFORM),mac_ifort)
  F90 = mpif90
  F77 = mpif77
  OPTS = -O3 -fpp -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'" -parallel -heap-arrays
  FOPTS = -O3 
  LIBS = -lmpi -llapack
  DEST = /Users/cjweiss/bin
  MANPATH = /usr/local/man/man1/.
endif

ifeq ($(PLATFORM),anodyne)
  F90 = mpif90
  F77 = mpif90
  OPTS = -O3 -fpp -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'" -openmp -parallel -heap-arrays
  OPTS = -O3 -fpp -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'" -heap-arrays
  FOPTS = -O3
  LIBS =  -lmpi 
  DEST = /Users/cjweiss/bin
  MANPATH = /usr/local/man/man1/.
endif

ifeq ($(PLATFORM),altix)
  F90 = ifort
  F77 = ifort
  OPTS =   -O2 -ipo -shared-intel -fpp -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'"
  FOPTS = -O2  -ipo
  LIBS = -lmpi -lmkl_lapack  -lmkl_ipf 
  DEST = /home/cjweiss/bin/.
  MANPATH = /usr/local/man/man1/.
endif

ifeq ($(PLATFORM),cabernet) 
  F90 = gfortran
  F77 = gfortran
  F90 = mpif90
  F77 = mpif77
  FOPTS =   -O3 
  OPTS =   -O3 -x f95-cpp-input -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'"
  INCL = -I/usr/include
  LIBS = -lmpi -llapack
  DEST = /home/cjweiss/bin/.
  MANPATH = /usr/local/man/man1/.
endif

ifeq ($(PLATFORM),hess) 
  F90 = ifort
  F77 = ifort
  F90 = mpif90
  F77 = mpif77
  OPTS = -O3 -i_dynamic -fpp -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'"
  LIBS = -lmpi -lmkl_lapack -lmkl 
  DEST = /home/cjweiss/bin
  MANPATH = /usr/local/man/man1/.
endif

ifeq ($(PLATFORM),mac_gfortran) 
  F90 = mpif90
  F77 = mpif90
  FOPTS =   -O3 
  OPTS =  -fcheck=bounds -m64 -O3 -x f95-cpp-input -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'" -fopenmp
  OPTS =  -m64 -O3 -x f95-cpp-input -DSVN_VERSION_DEF="svn_version_string = '`svnversion`'" -fopenmp -lmpi
  INCL = -I/usr/include
  LIBS = -llapack
  DEST = /home/cjweiss/bin/.
  MANPATH = /usr/local/man/man1/.
endif

#############################
#  Name for the executable  #
#############################

SERIAL_EXEC	= aphid_serial
PARALLEL_EXEC	= aphid_parallel

####################################
#  Destination for the executable  #
####################################

DEST		= ~/bin/.

#####################################
#  Object files for the executable  #
#####################################

mods		= aphid_mod.o
main		= aphid.o
qmr             = aphid_qmr.o

serial		= aphid_serial.o
parallel	= aphid_parallel.o

support		= $(qmr) $(main)

####################################
#  Build lines for the executable  #
####################################

serial:	 	$(mods) $(support) $(serial) 
	$(F90)  $(mods) $(support) $(serial) $(LIBS) $(OPTS) $(INCL) -o $(DEST)/$(SERIAL_EXEC)

parallel: 	$(mods) $(support) $(parallel)
	$(F90)  $(mods) $(support) $(parallel) $(LIBS) $(OPTS) $(INCL) -o $(DEST)/$(PARALLEL_EXEC)
	cp aphid_mpi $(DEST)/.

man:
	sudo cp ../man/aphid.1 $(MANPATH)/.

################################
#  Rules for Creating Objects  #
################################

.SUFFIXES: .f90 .o .h90

.f90.o:	
	$(F90)	$(OPTS) $(INCL) -c -o $@ $< 

.SUFFIXES: .f

.f.o:
	$(F77) $(OPTS) $(INCL) -c -o $@ $<

#################
#  Maintenance  #
#################

.PHONY:	clean
clean:
	rm -f *.o *.mod

