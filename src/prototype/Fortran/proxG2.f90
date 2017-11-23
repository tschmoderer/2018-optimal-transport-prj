!! etant donnée 4 variable mbar,fbar, m et f
!! calcul l'opérateur de proximité de G2
!! mbart, fbart, mt,ft

subroutine proxG2(U,V,Ut,Vt) 
    implicit none
	include 'global.inc'
	
	double precision, dimension((N+2)*(Q+1)+(N+1)*(Q+2)), intent(in) :: U;
	double precision, dimension((N+1)*(Q+1)+(N+1)*(Q+1)), intent(in) :: V;
	double precision, dimension((N+2)*(Q+1)+(N+1)*(Q+2)), intent(out) :: Ut;
	double precision, dimension((N+1)*(Q+1)+(N+1)*(Q+1)), intent(out) :: Vt;
	
	Ut = matmul(pG2,U+matmul(transpose(Interp),V));
	Vt = matmul(Interp,Ut);
end subroutine
