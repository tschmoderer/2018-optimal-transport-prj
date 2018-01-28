program test 
    implicit none 
    integer, parameter :: N = 200, Q = 300;
    double precision, dimension(1:Q+1,1:N+1,2) :: ASrArw, rw
    double precision, dimension(1:Q+3,1:N+1) :: Arw, rArw
    double precision, dimension(1:Q+1,1:N+1) :: rm, rdxrm, dxrm, dxSrdxrm
    double precision, dimension(1:Q+1,1:N+1) :: rf, rdtrf, dtrf, dtSrdtrf

    double precision, dimension(1:Q+3,1:N+1) :: resh1, resh2 
    double precision, dimension((Q+3)*(N+1)) :: flat1, flat2 
! Variables test projection
    double precision, dimension(1:Q+1,1:N+1,2) :: w
    double precision, dimension(1:Q+3,1:N+1) :: Aw 
    double precision, dimension(1:Q+3,1:N+1) :: y
    double precision :: epsilon = 1e-10
    double precision, dimension(1:N+1) :: f0,f1
    double precision, dimension(1:(Q+3)*(N+1)) :: b,r
    double precision, dimension(1:Q+1,1:N+1,2) :: pC
    double precision, dimension(1:Q+3,1:N+1) :: sol

    call random_number(rw)
    call random_number(rArw)
    call random_number(rm)
    call random_number(rdxrm)
    call random_number(rf)
    call random_number(rdtrf)

    call dx(dxrm,rm,N,Q)
    call dxS(dxSrdxrm,rdxrm,N,Q)

    call dt(dtrf,rf,N,Q)
    call dtS(dtSrdtrf,rdtrf,N,Q)

    print *, 'divergence m : ', sum(rm*dxSrdxrm - dxrm*rdxrm)
    print *, 'dérivée en t : ', sum(rf*dtSrdtrf - dtrf*rdtrf)

    call A(Arw,rw,N,Q);
    call AS(ASrArw,rArw,N,Q)

    print *, 'AAS   : ' , sum(rw*ASrArw) - sum(rArw*Arw)

    call random_number(resh1) 
    call flat(flat1,resh1,N,Q)
    call resh(resh2,flat1,N,Q)
    print *, 'Reshape test : ', sum(resh1 - resh2)

    call random_number(flat1)
    call resh(resh1, flat1,N,Q) 
    call flat(flat2,resh1,N,Q) 
    print *, ' Flat test :', sum(flat1 -flat2)

! Test projection 

    call random_number(w)

    call gauss(f0,0.5,0.05,N)
    call gauss(f1,0.8,0.05,N)

    call normalise(epsilon + f0,f0,N)
    call normalise(epsilon + f1,f1,N)

    y = 0;
    y(Q+2,:) = f0
    y(Q+3,:) = f1
    call A(Aw,w,N,Q)
!    print *, 'test A >0', sum(w*Aw)
    print *, 'Error avant projection :', sum(dabs(Aw -y))/sum(y)
return
    call flat(b,y-Aw,N,Q)
    call cg(r,b,N,Q)

    call resh(sol,r,N,Q)
    call AS(pC,sol,N,Q)

    pC = w + pC;

    call A(Aw,pC,N,Q)
    print *, 'Error post projection :', sum(dabs(Aw -y))/sum(y)
end program

subroutine cg(x,b,N,Q) 
    integer, intent(in) :: N,Q
    double precision, dimension(1:(Q+3)*(N+1)), intent(in) :: b
    double precision, dimension(1:(Q+3)*(N+1)), intent(out) :: x
    double precision, dimension(1:(Q+3)*(N+1)) :: x0, x1, r0, r1, p0, p1, Ap0
    double precision, dimension(1:Q+1,1:N+1,2) :: tmpA
    double precision, dimension(1:Q+3,1:N+1) :: tmpX, tmpAS, tmpR, tmpP

    double precision :: alpha,beta
    integer :: k, niter = 100

    x0 = 0;

    ! flat(A(AS(resh(x0))))
    call resh(tmpX,x0,N,Q)
    call AS(tmpA,tmpX,N,Q)
    call A(tmpR,tmpA,N,Q)
    call flat(r0,tmpR,N,Q)

    r0 = b - r0;
    p0 = r0;

    do while (k .LT. niter) 
      ! flat(A(AS(resh(p0))))
      call resh(tmpX,p0,N,Q)
      call AS(tmpA,tmpX,N,Q)
      call A(tmpR,tmpA,N,Q)
      call flat(Ap0,tmpR,N,Q)

      alpha = sum(r0*r0)/sum(p0*Ap0)
      r1 = r0 - alpha*Ap0
      x1 = x0 + alpha*p0

      beta = sum(r1*r1)/sum(r0*r0)
      p1 = r1 + beta*p0

      if (sum(dabs(r1)) .LT. 1e-15) then 
        exit; 
      end if
      print *, 'iteration CG ', k, sum(dabs(r1))
      r0 = r1;
      p0 = p1;
      x0 = x1;
      k = k+1;
    end do
    x = x1;
end subroutine

subroutine A(Aw,w,N,Q) 
    integer, intent(in) :: N,Q;
    double precision, dimension(1:Q+1,1:N+1,2), intent(in) :: w 
    double precision, dimension(1:Q+3,1:N+1), intent(out) :: Aw
    double precision, dimension(1:Q+1,1:N+1) :: tmpdiv, tmpdt

    call dxS(tmpdiv,w(:,:,1),N,Q)
    call dt(tmpdt,w(:,:,2),n,Q)

    Aw(1:Q+1,:) = -tmpdiv + tmpdt
    Aw(Q+2,:) = w(Q+1,:,2)
    Aw(Q+3,:) = w(1,:,2)
end subroutine

subroutine AS(ASAw,Aw,N,Q) 
    integer, intent(in) :: N,Q;
    double precision, dimension(1:Q+3,1:N+1), intent(in) :: Aw
    double precision, dimension(1:Q+1,1:N+1,2), intent(out) :: ASAw 
    double precision, dimension(1:Q+1,1:N+1) :: tmpdx, tmpdtS

    call dx(tmpdx,Aw(1:Q+1,:),N,Q)
    call dtS(tmpdtS,Aw(1:Q+1,:),N,Q)

    ASAw(:,:,1) = -tmpdx
    ASAw(:,:,2) = tmpdtS

    ASAw(1,:,2) = ASAw(1,:,2) + Aw(Q+3,:);
    ASAw(Q+1,:,2) = ASAw(Q+1,:,2) + Aw(Q+2,:);
end subroutine


subroutine dx(dm,m,N,Q) 
    integer, intent(in) :: N,Q;
    double precision, dimension(1:Q+1,1:N+1), intent(in) :: m;
    double precision, dimension(1:Q+1,1:N+1), intent(out) :: dm;
    
    dm = 0;
    dm(:,1:N) = N*(m(:,2:N+1) - m(:,1:N))
end subroutine

subroutine dt(df,f,N,Q) 
    integer, intent(in) :: N,Q;
    double precision, dimension(1:Q+1,1:N+1), intent(in) :: f;
    double precision, dimension(1:Q+1,1:N+1), intent(out) :: df;
    
    df = 0;
    df(1:Q,:) = Q*(f(2:Q+1,:) - f(1:Q,:))
end subroutine

subroutine dxS(m,dm,N,Q) 
    integer, intent(in) :: N,Q;
    double precision, dimension(1:Q+1,1:N+1), intent(in) :: dm;
    double precision, dimension(1:Q+1,1:N+1), intent(out) :: m;

    m(:,1)   = -dm(:,1)
    m(:,2:N) = dm(:,1:N-1) - dm(:,2:N)
    m(:,N+1) = dm(:,N)

    m(:,1:N+1) = N*m(:,1:N+1) 
end subroutine

subroutine dtS(f,df,N,Q) 
    integer, intent(in) :: N,Q;
    double precision, dimension(1:Q+1,1:N+1), intent(in) :: df;
    double precision, dimension(1:Q+1,1:N+1), intent(out) :: f;

    f(1,:)   = -df(1,:)
    f(2:Q,:) = df(1:Q-1,:) - df(2:Q,:)
    f(Q+1,:) = df(Q,:) 
    f(1:Q+1,:) = Q*f(1:Q+1,:)
end subroutine

subroutine flat(fx,x,N,Q) 
    integer, intent(in) :: N,Q;
    double precision, dimension(1:Q+3,1:N+1), intent(in) :: x
    double precision, dimension(1:(Q+3)*(N+1)), intent(out) :: fx
    fx = reshape(x,(/(Q+3)*(N+1)/)) ! une colonne pleins de lignes
end subroutine

subroutine resh(rx,x,N,Q)
    integer, intent(in) :: N,Q;
    double precision, dimension(1:(Q+3)*(N+1)), intent(in) :: x
    double precision, dimension(1:Q+3,1:N+1), intent(out) :: rx
    rx = reshape(x,(/Q+3,N+1/));
end subroutine

subroutine normalise(f,res,N)
    integer, intent(in) :: N
    double precision, dimension(1:N+1), intent(in) :: f
    double precision, dimension(1:N+1), intent(out) :: res
    res = f/sum(f)
end subroutine

subroutine gauss(f,mu,sigma,N)
    real, intent(in) :: mu, sigma
    integer, intent(in) :: N 
    double precision, dimension(1:N+1), intent(out) :: f
    integer :: i
    do i = 1,N+1 
        f(i) = exp(-0.5*(((i-1)/(1.0*N) - mu)/sigma)**2)
    end do
end subroutine