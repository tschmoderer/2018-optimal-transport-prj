subroutine proxG2(mbart,fbart,mt,ft,mbar,fbar,m,f,N,Q) 
    implicit none
    
    integer, intent(in) :: N,Q;
    real, dimension(N+1,Q+1), intent(in) :: m,f;
    real, dimension(N+2,Q+1), intent(in) :: mbar;
    real, dimension(N+1,Q+2), intent(in) :: fbar;
    real, dimension(N+1,Q+1), intent(out) :: mt,ft;
    real, dimension(N+2,Q+1), intent(out) :: mbart;
    real, dimension(N+1,Q+2), intent(out) :: fbart;

    real, dimension(N+1,N+2) :: Interpm;
    real, dimension(Q+2,Q+1) :: Interpf;
    real, dimension(N+2,N+1) :: InterpAdjM;
    real, dimension(Q+1,Q+2) :: InterpAdjF;
    real, dimension(N+1,N+1) :: Idm;
    real, dimension(N+1,N+1) :: Am;
    real, dimension(Q+1,Q+1) :: Idf;
    real, dimension(Q+1,Q+1) :: Af;
    real, dimension(N+2,Q+1) :: Bm;
    real, dimension(N+1,Q+2) :: Bf;

    integer :: i,j;
    integer :: IPIV,INFO;

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

    do i=1,N+1
        do j = 1,N+2
            if (i .EQ. j) then 
                Interpm(i,j) = 1;
            else if (i .EQ. j-1) then
                Interpm(i,j) = 1;
            end if
        end do 
    end do 
    Interpm = 0.5*Interpm;
    InterpAdjM = transpose(Interpm);
    
    do i=1,Q+2
        do j=1,Q+1
            if (i .EQ. j) then
                Interpf(i,j) = 1;
            else if (i .EQ. j+1) then 
                Interpf(i,j) = 1;
            end if
        end do
    end do
    Interpf = 0.5*Interpf;
    InterpAdjF = transpose(Interpf);

    do i=1,N+1
        do j=1,N+1
            if (i .EQ. j) then 
                Idm(i,j) = 1;
            end if 
        end do
    end do 

   do i=1,Q+1
        do j=1,Q+1
            if (i .EQ. j) then 
                Idf(i,j) = 1;
            end if 
        end do
    end do 


    !Am = Idm + matmul(Interpm,InterpAdjM);
    !Af = Idf + matmul(InterpAdjF,Interpf);

    !! second membre !! 

    !Bm = mbar + InterpAdjM*m;
    !Bf = fbar + f*InterpAdjF;

   ! call SGESV(N+1,N+1,Am,N+1,IPIV,Bm,N+1,INFO);
   ! call SGESV(N+1,N+1,Af,N+1,IPIV,Bf,N+1,INFO);

    !mbart = Bm;
    !fbart = Bf;
    !mt = matmul(Interpm,mbart);
    !ft = matmul(fbart,Interpf);
end subroutine