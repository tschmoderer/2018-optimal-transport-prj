!! opérateur de proximité G1
!! in : mbar, fbar, m, f
!! out : Umbar Ufbar, Vm,Vf

subroutine proxG1(U,V,pU,pV)
    implicit none
    include 'global.inc'

	double precision, dimension((N+2)*(Q+1)+(N+1)*(Q+2)), intent(in) :: U;
	double precision, dimension((N+1)*(Q+1)+(N+1)*(Q+1)), intent(in) :: V;
	double precision, dimension((N+2)*(Q+1)+(N+1)*(Q+2)), intent(out) :: pU;
	double precision, dimension((N+1)*(Q+1)+(N+1)*(Q+1)), intent(out) :: pV;
	
    call projC(U,pU);
    call proxJ(V,pV);
end subroutine
