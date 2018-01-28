FC = gfortran
FFLAGS  = -g -Wall
EXEC = transport

.DEFAULT_GOAL = all
.PHONY = clean

SRCS  =  $(wildcard ∗.f90)
OBJS = $(SRCS:.f90=.o)

%.o : %.f90
	$(FC) $(FFLAGS) −c $<

all : projC.o transport.o
	$(FC) $(FFLAGS) -o $(EXEC) *.o

clean : 
	rm *.out *.o $(EXEC) 