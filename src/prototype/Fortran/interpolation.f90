!! Etant donné des variables décentrées mbar,fbar !!
!! interpol les valeurs sur la grille centrée m,f !!

subroutine interpolation(Interp)
    implicit none
    include 'global.inc'
  
    integer :: i,j;

    double precision, dimension((N+1)*(Q+1),(N+2)*(Q+1)) :: Interpm;
    double precision, dimension((N+1)*(Q+1),(N+1)*(Q+2)) :: Interpf;
    double precision, dimension(Q+1,Q+2) :: tmp;

	Interpm = 0;
	Interpf = 0;
	tmp = 0;
	
	! Construction matrice d'interpolation pour m !
    do i=1,(N+1)*(Q+1) ! parcours les lignes
        do j = 1,(N+2)*(Q+1) ! parcours les colonnes 
            if (i .EQ. j) then 
                Interpm(i,j) = 0.5;
            else if (j .EQ. i + Q + 1) then
                Interpm(i,j) = 0.5;
            end if
        end do 
    end do 
    
    do i = 1,Q+1
		do j = 1,Q+2
			if (i .EQ. j) then
				tmp(i,j) = 1;
			else if (j .EQ. i+1) then
				tmp(i,j) = 1;
			end if
		end do
	end do

    ! construction matrice d'interpolation pour f ! 
    do i = 1,N+1
		Interpf((i-1)*(Q+1)+1:i*(Q+1),(i-1)*(Q+2)+1:i*(Q+2)) = tmp;
    end do
    
    Interp(1:(N+1)*(Q+1),1:(N+2)*(Q+1)) = Interpm;
    Interp((N+1)*(Q+1)+1:2*(N+1)*(Q+1),(N+2)*(Q+1)+1:2*(N+1)*(Q+1)+N+Q+2) = Interpf;

end subroutine
