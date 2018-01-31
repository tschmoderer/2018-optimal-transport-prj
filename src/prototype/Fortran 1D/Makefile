FC = gfortran
FFLAGS  = -g -Wall -Wno-tabs
EXEC = transport

.DEFAULT_GOAL = all
.PHONY = clean

SRCS  =  $(wildcard ∗.f90)
OBJS = $(SRCS:.f90=.o)

all : procedures.f90 transport.f90
	$(FC) $(FFLAGS) -o $(EXEC) procedures.f90 transport.f90

test : procedures.f90 test.f90
	$(FC) $(FFLAGS) procedures.f90 test.f90

new : procedures.f90 new.f90
	$(FC) $(FFLAGS) procedures.f90 new.f90 -o new

clean : 
	rm *.out *.o *.mod $(EXEC) 
