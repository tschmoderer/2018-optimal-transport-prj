!! etant donnée 4 variable mbar,fbar, m et f
!! calcul l'opérateur de proximité de G2
!! mbart, fbart, mt,ft

subroutine proxG2(mbart,fbart,mt,ft,mbar,fbar,m,f,N,Q) 
    implicit none
    
    integer, intent(in) :: N,Q;
    double precision, dimension(Q+1,N+1), intent(in) :: m,f;
    double precision, dimension(Q+1,N+2), intent(in) :: mbar;
    double precision, dimension(Q+2,N+1), intent(in) :: fbar;
    double precision, dimension(Q+1,N+1), intent(out) :: mt,ft;
    double precision, dimension(Q+1,N+2), intent(out) :: mbart;
    double precision, dimension(Q+2,N+1), intent(out) :: fbart;

    double precision, dimension(N+2,N+1) :: Interpm;
    double precision, dimension(Q+1,Q+2) :: Interpf;
    double precision, dimension(N+1,N+2) :: InterpAdjM;
    double precision, dimension(Q+2,Q+1) :: InterpAdjF;
    double precision, dimension(Q+2,Q+2) :: Idf;
    double precision, dimension(Q+2,Q+2) :: Af;
    double precision, dimension(Q+2,N+1) :: Bf;
    double precision, dimension(N+2,N+2) :: Idm;
    double precision, dimension(N+2,N+2) :: Am;
    double precision, dimension(N+2,Q+1) :: Bm;

    integer :: i,j;
    integer :: INFO;
    integer, dimension(N+2) :: IPIVm;
    integer, dimension(Q+2) :: IPIVf;

    external SGESV
  
    Interpm = 0;
    Interpf = 0;
    InterpAdjM = 0;
    InterpAdjF = 0;
    Idm = 0;
    Am = 0;
    Idf = 0;
    Af = 0;
    Bm = 0;
    Bf = 0;

	! Construction matrice d'interpolation pour m !
    do i=1,N+2
        do j = 1,N+1
            if (i .EQ. j) then 
                Interpm(i,j) = 0.5;
            else if (j .EQ. i-1) then
                Interpm(i,j) = 0.5;
            end if
        end do 
    end do 

    ! construction matrice d'interpolation pour f ! 
    do i=1,Q+1
        do j=1,Q+2
            if (i .EQ. j) then
                Interpf(i,j) = 0.5;
            else if (j .EQ. i+1) then 
                Interpf(i,j) = 0.5;
            end if
        end do
    end do
 
    InterpAdjM = transpose(Interpm); ! adjoint de Im
    InterpAdjF = transpose(Interpf); ! adjoint de If

    do i=1,N+2
        do j=1,N+2
            if (i .EQ. j) then 
                Idm(i,j) = 1;
            end if 
        end do
    end do 

   do i=1,Q+2
       do j=1,Q+2
            if (i .EQ. j) then 
                Idf(i,j) = 1;
            end if 
       end do
    end do 

    Am = transpose(Idm + matmul(Interpm,InterpAdjM));
    Af = Idf + matmul(InterpAdjF,Interpf);

    !! second membre !! 

   Bm = transpose(mbar + matmul(m,InterpAdjM));
   Bf = fbar + matmul(InterpAdjF,f);

   call SGESV(N+2,Q+1,Am,N+2,IPIVm,Bm,N+2,INFO);
   call SGESV(Q+2,N+1,Af,Q+2,IPIVf,Bf,Q+2,INFO);


    mbart = transpose(Bm);
    fbart = Bf;
    mt = matmul(mbart,Interpm);
    ft = matmul(Interpf,fbart);
end subroutine
