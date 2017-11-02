subroutine check_adjoint(N,Q)
    implicit none
    integer, intent(in) :: N,Q;

    !! Variables tests Boundary!!
    
    double precision, dimension(Q+1,N+2) :: mbar,boundAdjMres;
    double precision, dimension(Q+2,N+1) :: fbar,boundAdjFres;
    double precision, dimension(Q+1) :: mleft,mright,boundResLeft,boundResRight;
    double precision, dimension(N+1) :: fup,fdown,boundResUp,boundResDown;

    !! Variables pour l'interpolation !!
    
    double precision, dimension(Q+1,N+2) :: mbarInterpAdj;
    double precision, dimension(Q+2,N+1) :: fbarInterAdj;
    double precision, dimension(Q+1,N+1) :: m,f,mInterp,fInterp;

    !! Variables pour la divergence !! 

    double precision, dimension(Q+1,N+2) :: divAdjMres;
    double precision, dimension(Q+2,N+1) :: divAdjFres;
    double precision, dimension(Q+1,N+1) :: d,dRes;

    !! Variables pour tous le monde !!
    double precision :: s1,s2;


    !! Initialisation Boundary !!
    boundAdjMres = 0; boundAdjFres = 0;
    boundResLeft = 0; boundResRight = 0;
    boundResUp = 0; boundResDown = 0;
    call RANDOM_NUMBER(mbar);
    call RANDOM_NUMBER(fbar);
    call RANDOM_NUMBER(mleft);
    call RANDOM_NUMBER(mright);
    call RANDOM_NUMBER(fup);
    call RANDOM_NUMBER(fdown);

   
    print *, "Check Boundary adjoint";
    call boundary(mbar,fbar,boundResLeft,boundResRight,boundResUp,boundResDown,N,Q);
    call boundary_adjoint(boundAdjMres,boundAdjFres,mleft,mright,fup,fdown,N,Q);

    s1 = sum(boundResLeft*mleft) + sum(boundResRight*mright) + sum(boundResUp*fup) + sum(boundResDown*fdown);
    s2 = sum(boundAdjFres*fbar) + sum(boundAdjMres*mbar);

    print *, abs(s1-s2)/s1;

    !! Initialisation divergence !!
    divAdjMres = 0;
    divAdjFres = 0;
    d = 0;    
    call RANDOM_NUMBER(dRes);

    print *, "Check divergence adjoint";

    call divergence(mbar,fbar,d,N,Q);
    call divergence_adjoint(divAdjMres,divAdjFres,dRes,N,Q);

    s1 = sum(dRes*d);
    s2 = sum(divAdjMres*mbar) + sum(divAdjFres*fbar);

    print *, abs(s1-s2)/(1.0*s1);
    
 !! Initialisation Interpolation !!

    mbarInterpAdj = 0;
    fbarInterAdj = 0;
    mInterp = 0;
    fInterp = 0;
    call RANDOM_NUMBER(m);
    call RANDOM_NUMBER(f);

    print *, "Check interpolation adjoint";
  
    call interpolation(mbar,fbar,mInterp,fInterp,N,Q);
    call interpolation_adjoint(mbarInterpAdj,fbarInterAdj,m,f,N,Q);  

    s1 = sum(mInterp*m) + sum(fInterp*f);
    s2 = sum(mbarInterpAdj*mbar) + sum(fbarInterAdj*fbar);

    print *, abs(s1-s2)/(1.0*s1);

end subroutine