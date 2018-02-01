FC = gfortran
FFLAGS  = -g -Wall -Wno-tabs
EXEC = transport

.DEFAULT_GOAL = all
.PHONY = clean

SRCS  =  $(wildcard ∗.f90)
OBJS = $(SRCS:.f90=.o)

all : transport.f90
	$(FC) $(FFLAGS) -o $(EXEC) transport.f90


clean : 
	rm *.out *.o *.mod $(EXEC)
