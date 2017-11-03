subroutine DR(alpha,beta,gamma,N,Q)
    implicit none

    integer, intent(in) :: N,Q;
    double precision, intent(in) :: alpha,beta,gamma; 

    double precision, dimension(Q+1,N+1) :: zm,zf,wm0,wf0,wm1,wf1,wmtmp,wftmp;
    double precision, dimension(Q+1,N+2) :: zmbar,wmbar0,wmbar1,wmbartmp;
    double precision, dimension(Q+2,N+1) :: zfbar,wfbar0,wfbar1,wfbartmp;

    double precision, dimension(0:N) :: GcX;
    double precision, dimension(0:Q) :: GcT; 
    integer :: i,j,k;
    character(10) :: charI;
    character(200) :: cmd;
!    double precision :: R;

    !! initialisation
    do i=0,N
     GcX(i) = i/(1.0*N);
    end do
    do i =0,Q 
        GcT(i) = i/(1.0*Q);
    end do

    wm0 = 0; wm1 = 0;
    wf0 = 0; wf1 = 0;
    wmbar0 = 0; wmbar1 = 0;
    wfbar0 = 0; wfbar1 = 0;
!    call random_number(wmbar0);
!    call random_number(wfbar0);

!    call finitial(GcX,N,wfbar0(:,Q+2));
!    call ffinal(GcX,N,wfbar0(:,1));

    !! boucle principale
    do i = 0,50
        !! Svg des valeurs au début de la boucle 
        wmtmp = wm0; wftmp = wf0;
        wmbartmp = wmbar0; wfbartmp = wfbar0;

        !! z = ProxG2(w) !!
        call proxG2(zmbar,zfbar,zm,zf,wmbar0,wfbar0,wm0,wf0,N,Q);
        write(charI,'(I5.5)'), i
        open(7,file='results/f_'//trim(charI)//'.dat'); 
        write(7,*) "# ", "X ", "T ", "Z";
        do k = 1,Q+1 
            do j = 1,N+1 
                write(7,*) GcX(j-1), GcT(k-1), zf(k,j); 
            end do 
        end do
        close(7);
        write(cmd,*) 'set contour; set dgrid3d ', N, ',', Q, ';', 'splot "results/f_'//trim(charI), &
            '.dat" with lines; set title "Iteration Nb '//trim(charI),' "';
        open(8,file="plot.gnu"); write(8,*) cmd; close(8);
        call system("gnuplot plot.gnu");

        print *, "-3 : |wm1 - wm0| = ", sum(wm1-wm0);
        print *, "-2 : |wf1 - wf0| = ", sum(wf1-wf0);
        print *, "-1 : |wmbar1 - wmbar0| = ", sum(wmbar1-wmbar0);
        print *, "0 : |wfbar1 - wfbar0| = ", sum(wfbar1-wfbar0);
        
        !! RPrxG1 : 
    !    call proxG1(wmbar1,wfbar1,wm1,wf1,wmbar0,wfbar0,wm0,wf0,gamma,N,Q)
        call projC(wmbar1,wfbar1,wmbar0,wfbar0,N,Q);
        call proxJ(wm1,wf1,wm0,wf1,gamma,N,Q);  
     
     
        print *, "1 : |wm1 - wm0| = ", sum(wm1-wm0);
        print *, "2 : |wf1 - wf0| = ", sum(wf1-wf0);
        print *, "3 : |wmbar1 - wmbar0| = ", sum(wmbar1-wmbar0);
        print *, "4 : |wfbar1 - wfbar0| = ", sum(wfbar1-wfbar0);
       ! do while (1 .EQ. 1) 
       ! end do 
    !        print *, "Projection sur C mbar : ";
    !        do j = j,Q+1
    !            print *, wmbar1(:,j), "ENDL";
    !        end do 
    !        print *, "Projection sur C fbar : ";
    !        do j = 1,Q+2
    !            print *, wfbar1(:,j), "ENDL";
    !        end do
        !    do while (1 .EQ. 1)
        !    end do  
        
        wmbar0 = 2.0*wmbar1 - wmbar0; wfbar0 = 2.0*wfbar1 - wfbar0;
        wm0 = 2.0*wm1 - wm0; wf0 = 2.0*wf1 - wf0;

        !!RProxG2 : 
        call proxG2(wmbar1,wfbar1,wm1,wf1,wmbar0,wfbar0,wm0,wf0,N,Q);
        wmbar0 = 2.0*wmbar1 - wmbar0; wfbar0 = 2.0*wfbar1 - wfbar0;
        wm0 = 2.0*wm1 - wm0; wf0 = 2.0*wf1 - wf0;

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

!    open(0,file="solution/m"); write(0,*), zm; close(0);
!    open(1,file="solution/f"); write(1,*), zf; close(1);
!    open(2,file="solution/parameters"); write(2,*), N, Q; close(2);

    print *, "Fin de l'algorithme";
    stop;

!    print *, "zf : "
!    do i=1,N+1
!        print *, sum(zf(:,i)), "NN";
!    end do
!    print *, "zm : "
!    do i=1,N+1
!        print *, sum(zm(:,i)), "NN";
!    end do
!    print *, "zfbar : "
!    do i=1,Q+2
!        print *, sum(zfbar(:,i)), "NN";
!    end do
!    print *, "zmbar : "
!    do i=1,Q+1
!        print *, sum(zmbar(:,i)), "NN";
!    end do
    
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