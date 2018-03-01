# Copyright (C) 1997-2001 Imperial Cancer Research Technology
# author: Gidon Moont

# Biomolecular Modelling Laboratory
# Imperial Cancer Research Fund
# 44 Lincoln's Inn Fields
# London WC2A 3PX

# +44 (0)20 7269 3348
# http://www.bmm.icnet.uk/

#############

FFTW_DIR        = /path/to/FFTW/2.1.5/

#############

# You may need/want to edit some of these
#
# Hint: For the CC_FLAGS have a look at what the fftw build used

SHELL           = /bin/sh

# Use of the OpenMPI compiler, change it according to your MPI environment
CC		= mpicc

# Something more aggresive than -03 produces numerical unstability
CC_FLAGS	= -DUSE_MPI -O3 -pthread

CC_LINKERS      = -lm

STRIP           = strip

SECURITY	= chmod 555

####################################################

# You should not be editing anything below here

CC_FLAGS_FULL   = -I$(FFTW_DIR)/include $(CC_FLAGS)

# Using double precision here
FFTW_LINKERS    = -L$(FFTW_DIR)/lib -ldrfftw -ldfftw

# Old configuration, might be useful in some old machines
#CC_FLAGS_FULL	= -I$(FFTW_DIR)/fftw -I$(FFTW_DIR)/rfftw $(CC_FLAGS)
#FFTW_LINKERS    = -L$(FFTW_DIR)/fftw/.libs -L$(FFTW_DIR)/rfftw/.libs -lrfftw -lfftw

#############

.SUFFIXES:	.c .o

.c.o:
		$(CC) $(CC_FLAGS_FULL) -c $<

#############

LIBRARY_OBJECTS = manipulate_structures.o angles.o coordinates.o electrostatics.o grid.o qsort_scores.o

PROGRAMS = ftdock rpscore rpdock filter build centres randomspin

all:		$(PROGRAMS)
 
#############

ftdock:		ftdock.o $(LIBRARY_OBJECTS) structures.h
		$(CC) $(CC_FLAGS_FULL) -o $@ ftdock.o $(LIBRARY_OBJECTS) $(FFTW_LINKERS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@

#############

rpdock:		rpdock.o $(LIBRARY_OBJECTS) structures.h
		$(CC) $(CC_FLAGS) -o $@ rpdock.o $(LIBRARY_OBJECTS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@


#############

rpscore:	rpscore.o $(LIBRARY_OBJECTS) structures.h
		$(CC) $(CC_FLAGS) -o $@ rpscore.o $(LIBRARY_OBJECTS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@

#############

filter:		filter.o $(LIBRARY_OBJECTS) structures.h
		$(CC) $(CC_FLAGS) -o $@ filter.o $(LIBRARY_OBJECTS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@

#############

build:		build.o $(LIBRARY_OBJECTS) structures.h
		$(CC) $(CC_FLAGS) -o $@ build.o $(LIBRARY_OBJECTS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@

#############

centres:	centres.o $(LIBRARY_OBJECTS) structures.h
		$(CC) $(CC_FLAGS) -o $@ centres.o $(LIBRARY_OBJECTS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@

#############

randomspin:	randomspin.o $(LIBRARY_OBJECTS) structures.h
		$(CC) $(CC_FLAGS) -o $@ randomspin.o $(LIBRARY_OBJECTS) $(CC_LINKERS)
		$(STRIP) $@
		$(SECURITY) $@

#############


clean:
		rm -f *.o core $(PROGRAMS)

#############

# dependencies

ftdock.o:			structures.h
rpscore.o:			structures.h
rpdock.o:			structures.h
filter.o:			structures.h
build.o:			structures.h
centres.o:			structures.h
randomspin.o:			structures.h

angles.o:			structures.h
coordinates.o:			structures.h
electrostatics.o:		structures.h
grid.o:				structures.h
manipulate_structures.o:	structures.h
qsort_scores.o:			structures.h
