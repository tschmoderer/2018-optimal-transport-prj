FC = gfortran
FFLAGS  = -g -o
LLFLAGS  = -llapack
EXEC = transport

.DEFAULT_GOAL = all
.PHONY = clean

SRCS  =  $(wildcard ∗.f90)
OBJS = $(SRCS:.f90=.o)

%.o : %.f90
	$(FC) $(FFLAGS) −c $<

all : *.f90	
	$(FC) $(FFLAGS) $(EXEC) *.f90 $(LLFLAGS)

clean : 
	rm *.out *.o $(EXEC) 		