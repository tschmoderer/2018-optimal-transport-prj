program new
    implicit none
    integer, parameter :: N = 10, Q = 10, niter = 200
    double precision, parameter :: eps = 1e-10, alpha = 1.0, g = 1.0
    double precision, dimension(N+1) :: f0, f1
    double precision, dimension(Q+1,N+1,2) :: z = 0, w0 = 0, w1 = 0
    double precision, dimension(niter) :: cout, minF, div
    integer :: i

    f0 = normalise(eps + gauss(0.5d0,0.05d0))
    f1 = normalise(eps + gauss(0.5d0,0.05d0))

    do i = 1,niter
        w1 = w0 + alpha*(proxJ(2*z-w0) - z)
        z = projC(w1) 

        cout(i) = J(z)
        print *, i, cout(i)
        minF(i) = minval(z(:,:,2))
    end do 

    contains

!! Gauss
    function gauss(mu,sigma) result(f) 
        double precision :: mu, sigma 
        double precision, dimension(N+1) :: f
        integer :: i
        do i = 1,N+1
            f(i) = exp(-0.5*((((i-1)/(1.0*N))-mu)/sigma)**2)
        end do 
    end function gauss

!! Normalise
    function normalise(f) result(nf) 
        double precision, dimension(N+1) :: f, nf
        nf = f/sum(f)
    end function

!! Cout 
    function J(w) result(c) 
        double precision, dimension(Q+1,N+1,2) :: w
        double precision :: c
        c = sum(w(:,:,1)**2/w(:,:,2))
    end function J 

!! Proximal de J 
    function proxJ(w) result(pw) 
        double precision, dimension(Q+1,N+1,2) :: w, pw
        double precision, dimension(Q+1,N+1) :: mt, ft, x0, x1, poly, dpoly
        integer :: k,l
        mt = w(:,:,1); ft = w(:,:,2);
        x0 = 1; x1 = 2; k = 0;

        do while (maxval(dabs(x0-x1)) .GT. 1e-5  .AND. k .LT. 1500)
            x0 = x1
            poly  = (x0-ft)*(x0+g)**2 - 0.5*g*mt**2
            dpoly = 2*(x0+g)*(x0-ft) + (x0+g)**2
            x1 = x0 - poly/dpoly
            k = k+1
        end do
        
        pw(:,:,2) = x1
        pw(:,:,1) = x1*mt/(x1+g) 

        do k = 1,Q+1
            do l = 1,N+1
                if (x1(k,l) .LT. 0) then 
                    pw(k,l,2) = eps
                    pw(k,l,1) = 0
                end if 
            end do
        end do 
    end function proxJ

!! Projection sur C 
    function projC(w) result(pc)
        double precision, dimension(Q+1,N+1,2) :: w, pc
        pc = 0
    end function projC
end program new