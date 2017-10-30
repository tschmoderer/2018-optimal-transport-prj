!! Calcul le cout d'une solution m, f sur al grille centrée

subroutine cost(R,m,f,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    double precision, dimension(N+1,Q+1), intent(in) :: m,f; 
    double precision, intent(out) :: R;

    integer :: i,j; 
    real :: x = 0.0;

    R = 0;
    do i = 1,N+1
        do j = 1,Q+1
            if (f(i,j) .GT. 0) then 
                R = R + m(i,j)**2/(1.0*f(i,j));
            else if ((f(i,j) .EQ. 0 ) .AND. (m(i,j) .EQ. 0)) then 
                R = R;
            else 
                R = R - log(x);
            end if
        end do 
    end do 
    
   !! second possibility more consicise but harder to handle all cases
   !R = 0;
   !R = sum(m**2/(1.0*f));
end subroutine