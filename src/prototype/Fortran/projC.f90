!! etant donné mbar et fbar, 
!! calcul la projection sur l'ensemble C 
!! Appel la résolution du pbm de poisson

subroutine projC(U,pU)
    implicit none
	include 'global.inc'
	double precision, dimension((N+2)*(Q+1)+(N+1)*(Q+2)), intent(in) :: U;
	double precision, dimension((N+2)*(Q+1)+(N+1)*(Q+2)), intent(out) :: pU;
	
	pU = matmul(P,U) + Cst;
	
end subroutine
