!! Etant donnée des valeurs sur la grille centrée m,f !!
!! Construit les valeurs décentrées mbar et fbar !!

subroutine interpolation_adjoint(mbar,fbar,m,f,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    
    real, dimension(N+2,Q+1), intent(out) :: mbar;
    real, dimension(N+1,Q+2), intent(out) :: fbar;
    real, dimension(N+1,Q+1), intent(in) :: m,f;
    real, dimension(N+1,N+2) :: Interpm;
    real, dimension(Q+2,Q+1) :: Interpf;
    real, dimension(N+2,N+1) :: InterpAdjM;
    real, dimension(Q+1,Q+2) :: InterpAdjF;

    integer :: i,j;
    !! test
    real, dimension(N+1,N+1) :: test1;
    real, dimension(Q+1,Q+1) :: test2;

    Interpm = 0;
    Interpf = 0;
    InterpAdjM = 0;
    InterpAdjF = 0;
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
    ! construction matrice d'interpolation de f
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

    mbar = matmul(InterpAdjM,m);
    fbar = matmul(f,InterpAdjF);
end subroutine