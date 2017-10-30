!! etant donnée 4 variable mbar,fbar, m et f
!! calcul l'opérateur de proximité de G2
!! mbart, fbart, mt,ft

subroutine proxG2(mbart,fbart,mt,ft,mbar,fbar,m,f,N,Q) 
    implicit none
    
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1), intent(in) :: m,f;
    double precision, dimension(N+2,Q+1), intent(in) :: mbar;
    double precision, dimension(N+1,Q+2), intent(in) :: fbar;
    double precision, dimension(N+1,Q+1), intent(out) :: mt,ft;
    double precision, dimension(N+2,Q+1), intent(out) :: mbart;
    double precision, dimension(N+1,Q+2), intent(out) :: fbart;

    double precision, dimension(N+1,N+2) :: Interpm;
    double precision, dimension(Q+2,Q+1) :: Interpf;
    double precision, dimension(N+2,N+1) :: InterpAdjM;
    double precision, dimension(Q+1,Q+2) :: InterpAdjF;
    double precision, dimension(N+2,N+2) :: Idm;
    double precision, dimension(N+2,N+2) :: Am;
    double precision, dimension(Q+2,Q+2) :: Idf;
    double precision, dimension(Q+2,Q+2) :: Af;
    double precision, dimension(N+2,Q+1) :: Bm;
    double precision, dimension(N+1,Q+2) :: Bf;

    integer :: i,j;
    integer :: INFO;
   ! integer, parameter :: LWORK = Q+2;
    double precision, dimension(Q+2) :: WORK;
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

    ! construction matrice d'interpolation de m
    do j=1,N+2
        do i = 1,N+1
            if (i .EQ. j) then 
                Interpm(i,j) = 0.5;
            else if (i .EQ. j-1) then
                Interpm(i,j) = 0.5;
            end if
        end do 
    end do 
    InterpAdjM = transpose(Interpm); ! son adjoint
    
  !   construction matrice d'interpolation de f
    do i=1,Q+2
        do j=1,Q+1
            if (i .EQ. j) then
                Interpf(i,j) = 0.5;
            else if (i .EQ. j+1) then 
                Interpf(i,j) = 0.5;
            end if
        end do
    end do
    InterpAdjF = transpose(Interpf); ! son adjointe

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

!    print *, "Interpf : ";
 !   do i=1,Q+1
  !  print *, Interpf(:,i), "NN";
   ! end do 

   ! print *, "InterpfAdj : ";
  !  do i =1,Q+2
 !   print *, InterpAdjF(:,i), "NN";
!    end do 

    Af = Idf + matmul(Interpf,InterpAdjF);
    !call DGETRI(Q+2,Af,Q+2,IPIVf,WORK,Q+2,INFO);
  Am = Idm + matmul(InterpAdjM,Interpm);
 !  Af = Idf + matmul(InterpAdjF,Interpf);

    !! second membre !! 

  Bm = mbar + matmul(InterpAdjM,m);
  Bf = fbar + matmul(f,InterpAdjF);

   ! print *, "Af : ";
  !  do i=1,Q+2
 !   print *, Af(:,i), "NN";
!    end do 

!    print *, "Bf : ";
 !   do i =1,Q+2
  !  print *, Bf(:,i), "NN";
   ! end do 
  
   call SGESV(N+2,Q+1,Am,N+2,IPIVm,Bm,N+2,INFO);
   call SGESV(Q+2,N+1,Af,Q+2,IPIVf,Bf,Q+2,INFO);

 !print *, "la solution Bf : ";
  !  do i =1,Q+2
   ! print *, Bf(:,i), "NN";
    !end do 
    mbart = Bm;
    fbart = Bf;
    mt = matmul(Interpm,mbart);
    ft = matmul(fbart,Interpf);
end subroutine