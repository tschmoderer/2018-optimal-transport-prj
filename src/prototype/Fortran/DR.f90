subroutine DR(alpha,beta,gamma,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    real, intent(in) :: alpha,beta,gamma; 

    double precision, dimension(N+1,Q+1) :: zm,zf,wm0,wf0,wm1,wf1,wmtmp,wftmp;
    double precision, dimension(N+2,Q+1) :: zmbar,wmbar0,wmbar1,wmbartmp;
    double precision, dimension(N+1,Q+2) :: zfbar,wfbar0,wfbar1,wfbartmp;

    double precision, dimension(0:N) :: GcX; ! test!
    integer :: i;

    !! initialisation
    do i=0,N
    GcX(i) = i/(1.0*N);
    end do
    wm0 = 1000; wm1 = 0;
    wf0 = 1000; wf1 = 0;
    wmbar0 = 1000; wmbar1 = 0;
    wfbar0 = 1000; wfbar1 = 0;
    !call finitial(GcX,N,wfbar0(:,Q+2));
   ! call random_number(wmbar0);

    do i = 0,100
        !! z = ProxG2(w) !!
        call proxG2(zmbar,zfbar,zm,zf,wmbar0,wfbar0,wm0,wf0,N,Q);
        !! w(l+1) = 
        wmtmp = wm0; wftmp = wf0;
        wmbartmp = wmbar0; wfbartmp = wfbar0;
        !! RPrxG1 : 
        call proxG1(wmbar1,wfbar1,wm1,wf1,wmbar0,wfbar0,wm0,wf0,gamma,N,Q)
    !    print *, "1 : |wm1 - wm0| = ", sum(wm1-wm0);
    !    print *, "2 : |wf1 - wf0| = ", sum(wf1-wf0);
    !    print *, "3 : |wmbar1 - wmbar0| = ", sum(wmbar1-wmbar0);
    !    print *, "4 : |wfbar1 - wfbar0| = ", sum(wfbar1-wfbar0);
        
        wmbar0 = 2*wmbar1 - wmbar0; wfbar0 = 2*wfbar1 - wfbar0;
        wm0 = 2*wm1 - wm0; wf0 = 2*wf1 - wf0;

        !!RProxG2 : 
        call proxG2(wmbar1,wfbar1,wm1,wf1,wmbar0,wfbar0,wm0,wf0,N,Q);
        wmbar0 = 2*wmbar1 - wmbar0; wfbar0 = 2*wfbar1 - wfbar0;
        wm0 = 2*wm1 - wm0; wf0 = 2*wf1 - wf0;

    !    print *, "5 : |wm1 - wm0| = ", sum(wm1-wm0);
    !    print *, "6 : |wf1 - wf0| = ", sum(wf1-wf0);
    !    print *, "7 : |wmbar1 - wmbar0| = ", sum(wmbar1-wmbar0);
    !    print *, "8 : |wfbar1 - wfbar0| = ", sum(wfbar1-wfbar0);


        wmbar1 = (1-0.5*alpha)*wmbartmp + 0.5*alpha*wmbar0;
        wfbar1 = (1-0.5*alpha)*wfbartmp + 0.5*alpha*wfbar0;
        wm1 = (1-0.5*alpha)*wmtmp + 0.5*alpha*wm0;
        wf1 = (1-0.5*alpha)*wftmp + 0.5*alpha*wf0;

        wm0 = wm1;wf0 = wf1;
        wmbar0 = wmbar1; wfbar0 = wfbar1;
    end do

!    print *, "zf : "
!    do i=1,N+1
!        print *, sum(zf(:,i)), "NN";
!    end do
    print *, "zm : "
    do i=1,N+1
        print *, sum(zm(:,i)), "NN";
    end do
!    print *, "zfbar : "
!    do i=1,Q+2
!        print *, sum(zfbar(:,i)), "NN";
!    end do
!    print *, "zmbar : "
!    do i=1,Q+1
!        print *, sum(zmbar(:,i)), "NN";
!    end do

    open(0,file="solution/m"); write(0,*), zm; close(0);
    open(1,file="solution/f"); write(1,*), zf; close(1);
    open(2,file="solution/parameters"); write(2,*), N, Q; close(2);
!    print *, "wf : "
!    do i=1,N+1
!        print *, sum(wf0(:,i)), "NN";
!    end do
!    print *, "wm : "
!    do i=1,N+1
!        print *, sum(wm0(:,i)), "NN";
!    end do
!    print *, "wfbar : "
!    do i=1,Q+2
!        print *, sum(wfbar0(:,i)), "NN";
!    end do
!    print *, "wmbar : "
!    do i=1,Q+1
!        print *, sum(wmbar0(:,i)), "NN";
!    end do
end subroutine