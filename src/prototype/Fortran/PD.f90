subroutine PD(alpha,beta,gamma,N,Q)
    implicit none

    integer, intent(in) :: N,Q;
    double precision, intent(in) :: alpha,beta,gamma; 
    ! variables d'itération
    double precision, dimension(Q+1,N+1) :: Vm0,Vm1,Vf0,Vf1;
    double precision, dimension(Q+1,N+2) :: Umbar0,Umbar1,Gmbar0,Gmbar1;
    double precision, dimension(Q+2,N+1) :: Ufbar0,Ufbar1,Gfbar0,Gfbar1;
    
    ! variables tampon 
    double precision, dimension(Q+1,N+1) :: V0m,V0f,tmp0,tmp1;
    double precision, dimension(Q+1,N+2) :: U0mbar, tmp2, tmp3;
    double precision, dimension(Q+2,N+1) :: U0fbar, tmp4, tmp5;

    double precision :: theta, sigma,tau;
    integer :: i; 

    !! Iniialisation : 
    Vm0 = 0; Vf0 = 1; Vm1 = 0; Vf1 = 0; 
    Umbar0 = 0; Ufbar0 = 0; Umbar1 = 0; Ufbar1 = 0;
    Gmbar0 = 0; Gfbar0 = 0; Gfbar1 = 0; Gfbar1 = 0;

    theta = 0.5; sigma = 85.0;
    tau = 0.99/(1.0*sigma*4);
    do i = 0,5
        !! Itération 1 : Prox G2*
        call interpolation(Gmbar0,Gfbar0,V0m,V0f,N,Q);
        Vm0 = Vm0 + sigma*Vm0; Vf0 = Vf0 + sigma*Vf0;
        call proxG2(tmp2,tmp4,Vm1,Vf1,tmp3,tmp5,Vm0/(1.0*sigma),Vf0/(1.0*sigma),N,Q);
        Vm1 = Vm0 - sigma*Vm1; Vf1 = Vf0 - sigma*Vf1;

        !! Itération 2: 
        call interpolation_adjoint(U0mbar,U0fbar,Vm1,Vf1,N,Q);
        Umbar0 = Umbar0 - tau*U0mbar;
        Ufbar0 = Ufbar0 - tau*U0fbar;
        call proxG1(Umbar1,Ufbar1,V0m,V0f,Umbar0,Ufbar0,tmp0,tmp1,tau,N,Q)

        !! Itération 3
        Gmbar1 = Umbar1 + theta*(Umbar1 - Umbar0);
        Gfbar1 = Ufbar1 + theta*(Ufbar1 - Ufbar0);

        !! MaJ des variables 
        Vm0 = Vm1; Vf0 = Vf1;
        Umbar0 = Umbar1; Ufbar0 = Ufbar1; 
        Gmbar0 = Gmbar1; Gmbar0 = Gmbar1;
    end do

end subroutine