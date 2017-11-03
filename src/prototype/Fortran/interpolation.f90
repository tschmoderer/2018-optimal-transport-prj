!! Etant donné des variables décentrées mbar,fbar !!
!! interpol les valeurs sur la grille centrée m,f !!

subroutine interpolation(mbar,fbar,m,f,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    integer :: i,j;
    double precision, dimension(Q+1,N+2), intent(in) :: mbar;
    double precision, dimension(Q+2,N+1), intent(in) :: fbar;
    double precision, dimension(Q+1,N+1), intent(out) :: m,f;
    double precision, dimension(N+2,N+1) :: Interpm;
    double precision, dimension(Q+1,Q+2) :: Interpf;

	Interpm = 0;
	Interpf = 0;
	! Construction matrice d'interpolation pour m !
    do i=1,N+2 ! parcours les lignes
        do j = 1,N+1 ! parcours les colonnes 
            if (i .EQ. j) then 
                Interpm(i,j) = 0.5;
            else if (j .EQ. i-1) then
                Interpm(i,j) = 0.5;
            end if
        end do 
    end do 

    ! construction matrice d'interpolation pour f ! 
    do i=1,Q+1 ! parcours lignes 
        do j=1,Q+2 ! parcours colonnes 
            if (i .EQ. j) then
                Interpf(i,j) = 0.5;
            else if (j .EQ. i+1) then 
                Interpf(i,j) = 0.5;
            end if
        end do
    end do

    m = matmul(mbar,Interpm);
    f = matmul(Interpf,fbar);
end subroutine
