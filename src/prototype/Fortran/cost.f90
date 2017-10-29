!! Calcul le cout d'une solution m, f sur al grille centrée

subroutine cost(R,m,f,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    real, dimension(N+1,Q+1), intent(in) :: m,f; 
    real, intent(out) :: R;

    R = sum(m**2/(1.0*f));
end subroutine