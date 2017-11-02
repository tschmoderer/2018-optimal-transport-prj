!! Etant donnée des valeurs sur la grille centrée m,f !!
!! Construit les valeurs décentrées mbar et fbar !!

subroutine interpolation_adjoint(mbar,fbar,m,f,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    
    double precision, dimension(Q+1,N+2), intent(out) :: mbar;
    double precision, dimension(Q+2,N+1), intent(out) :: fbar;
    double precision, dimension(Q+1,N+1), intent(in) :: m,f;
    double precision, dimension(N+2,N+1) :: Interpm;
    double precision, dimension(Q+1,Q+2) :: Interpf;
    double precision, dimension(N+1,N+2) :: InterpAdjM;
    double precision, dimension(Q+2,Q+1) :: InterpAdjF;

    integer :: i,j;
    
	Interpm = 0;
	Interpf = 0;

	! Construction matrice d'interpolation pour m !
    do i=1,N+2
        do j = 1,N+1
            if (i .EQ. j) then 
                Interpm(i,j) = 0.5;
            else if (j .EQ. i+1) then
                Interpm(i,j) = 0.5;
            end if
        end do 
    end do 

    ! construction matrice d'interpolation pour f ! 
    do i=1,Q+1
        do j=1,Q+2
            if (i .EQ. j) then
                Interpf(i,j) = 0.5;
            else if (i .EQ. j+1) then 
                Interpf(i,j) = 0.5;
            end if
        end do
    end do
 
    InterpAdjM = transpose(Interpm); ! adjoint de Im
    InterpAdjF = transpose(Interpf); ! adjoint de If

    mbar = matmul(m,InterpAdjM);
    fbar = matmul(InterpAdjF,f);
end subroutine