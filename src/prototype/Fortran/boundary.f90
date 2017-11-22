!! Construit la matrice de l'opérateur B !!
!! Timothée Schmoderer !! 
!! cc 2017!!

! function boundary
! implicit none
! include 'global.inc'
! double precision, dimension(size(B) :: boundary;

!boundary = B;
! end function
subroutine boundary(B)
    implicit none
    include 'global.inc'

    integer :: i,j;
    double precision, dimension(2*(Q+1),(N+2)*(Q+1)) :: Bm;
    double precision, dimension(2*(N+1),(Q+2)*(N+1)) :: Bf;  
    double precision, dimension(2,Q+2) :: tmp;

    ! Init ! 
    Bm = 0; Bf = 0; tmp = 0; B = 0;

    do i = 1,Q+1
        Bm(i,i) = 1;
        Bm(2*(Q+1)+1-i,(N+2)*(Q+1)+1-i) = 1;
    end do 

    tmp(1,1) = 1;
    tmp(2,Q+2) = 1;

    do i = 1,N+1
        Bf(2*i-1:2*i,(i-1)*(Q+2)+1:i*(Q+2)) = tmp;
    end do

    B(1:2*(Q+1),1:(N+2)*(Q+1)) = Bm;
    B(2*Q+3:2*Q+2*N+4,(N+2)*(Q+1)+1:2*(Q+1)*(N+1)+Q+N+2) = Bf;
end subroutine