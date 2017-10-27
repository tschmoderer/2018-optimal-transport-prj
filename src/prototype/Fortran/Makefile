FC = gfortran
FFLAGS  = -o
LDFLAGS  =
EXEC = transport

.DEFAULT_GOAL = all

SRCS  =  $(wildcard ∗.f90)
OBJS = $(SRCS:.f90=.o)

%.o : %.f90
	$(FC) $(FFLAGS) −c $<

all : transport.f90 finitial.f90 ffinal.f90 interpolation.f90 interpolation_adjoint.f90
	$(FC) $(FFLAGS) $(EXEC) transport.f90 finitial.f90 ffinal.f90 interpolation.f90 interpolation_adjoint.f90

clean : 
	rm *.out *.o $(EXEC) 