!! etant donné un champ de valeurs sur al grille centrée
!! construit des valeurs décentrées <=> -Grad

subroutine divergence_adjoint(mbar,fbar,d,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    real, dimension(N+2,Q+1), intent(out) :: mbar;
    real, dimension(N+1,Q+2), intent(out) :: fbar;
    real, dimension(N+1,Q+1), intent(in) :: d; 

    mbar(1,:) = -d(1,:);
    mbar(N+2,:) = d(N+1,:);
    mbar(2:N+1,:) = d(1:N,:) - d(2:N+1,:);

    fbar(:,1) = -d(:,1);
    fbar(:,Q+2) = d(:,Q+1);
    fbar(:,2:Q+1) = d(:,1:Q) - d(:,2:Q+1);
end subroutine