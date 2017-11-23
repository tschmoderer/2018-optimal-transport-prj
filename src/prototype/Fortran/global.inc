! Define all globals variable

integer, parameter :: N = 1 , Q = 1;
double precision, dimension(2*Q+2*N+4,2*(Q+1)*(N+1)+Q+N+2) :: B;
double precision, dimension((N+1)*(Q+1),2*(N+1)*(Q+1)+N+Q+2) :: D;
double precision, dimension(2*(N+1)*(Q+1),2*(N+1)*(Q+1)+N+Q+2) :: Interp;

double precision, dimension(2*(N+1)*(Q+1)+N+Q+2,2*(N+1)*(Q+1)+N+Q+2) :: P;
double precision, dimension(2*(N+1)*(Q+1)+N+Q+2) :: Cst;
double precision, dimension(2*(N+1)*(Q+1)+N+Q+2,2*(N+1)*(Q+1)+N+Q+2) :: pG2;

double precision :: g;
