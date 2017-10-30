subroutine DR(alpha,gamma,N,Q)
    implicit none
    integer, intent(in) :: N,Q;
    real, intent(in) :: alpha,gamma; 

    double precision, dimension(N+1,Q+1) :: zm,zf,wm0,wf0,wm1,wf1,wmtmp,wftmp;
    double precision, dimension(N+2,Q+1) :: zmbar,wmbar0,wmbar1,wmbartmp;
    double precision, dimension(N+1,Q+2) :: zfbar,wfbar0,wfbar1,wfbartmp;

    integer :: i;

    !! initialisation
    wm0 = 0; wm1 = 0;
    wf0 = 0; wf1 = 0;
    wmbar0 = 0; wmbar1 = 0;
    wfbar0 = 0; wfbar1 = 0;

    
    do i = 0,10
        !! z = ProxG2(w) !!
        call proxG2(zmbar,zfbar,zm,zf,wmbar0,wfbar0,wm0,wf0,N,Q);
        !! w(l+1) = 
        wmtmp = wm0; wftmp = wf0;
        wmbartmp = wmbar0; wfbartmp = wfbar0;
        !! RPrxG1 : 
        call proxG1(wmbar1,wfbar1,wm1,wf1,wmbar0,wfbar0,wm0,wf0,gamma,N,Q)
        wmbar0 = 2*wmbar1 - wmbar0; wfbar0 = 2*wfbar1 - wfbar0;
        wm0 = 2*wm1 - wm0; wf0 = 2*wf1 - wf0;
        !!RProxG2 : 
        call proxG2(wmbar1,wfbar1,wm1,wf1,wmbar0,wfbar0,wm0,wf0,N,Q);
        wmbar0 = 2*wmbar1 - wmbar0; wfbar0 = 2*wfbar1 - wfbar0;
        wm0 = 2*wm1 - wm0; wf0 = 2*wf1 - wf0;

        wmbar1 = (1-0.5*alpha)*wmbartmp + 0.5*alpha*wmbar0;
        wfbar1 = (1-0.5*alpha)*wfbartmp + 0.5*alpha*wfbar0;
        wm1 = (1-0.5*alpha)*wmtmp + 0.5*alpha*wm0;
        wf1 = (1-0.5*alpha)*wftmp + 0.5*alpha*wf0;

        wm0 = wm1;wf0 = wf1;
        wmbar0 = wmbar1; wfbar0 = wfbar1;
    end do
 call proxG2(zmbar,zfbar,zm,zf,wmbar0,wfbar0,wm0,wf0,N,Q);
    print *, "zf : "
    do i=1,N+1
        print *, sum(zf(:,i)), "NN";
    end do
    print *, "zm : "
    do i=1,N+1
        print *, sum(zm(:,i)), "NN";
    end do
    print *, "zfbar : "
    do i=1,N+1
        print *, sum(zfbar(:,i)), "NN";
    end do
    print *, "zmbar : "
    do i=1,N+1
        print *, sum(zmbar(:,i)), "NN";
    end do
end subroutine