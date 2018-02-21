clc 
clear all 


N = 9;
tmp = [1:N+1]; 
fprintf('Test 1D DCT : %f\n',sum(abs(dct(tmp) - mydct(tmp,N)))) 

dtmp = dct(tmp);
fprintf('Test 1D IDCT : %f\n',sum(abs(idct(dtmp) - myidct(dtmp,N)))) 


tmp = zeros(2,3);
tmp(1,1) = 1; tmp(1,2) = 2; tmp(1,3) = 3;
tmp(2,1) = 1; tmp(2,2) = 2; tmp(2,3) = 3;
tmp(3,1) = 1; tmp(3,2) = 2; tmp(3,3) = 3;
tmp(4,1) = 1; tmp(4,2) = 2; tmp(4,3) = 3;

fprintf('Test 2D DCT : %f\n',sum(sum(dct2(tmp) - mydct2(tmp,2,3)))) 

dtmp = dct2(tmp);

fprintf('Test 2D IDCT : %f\n',sum(sum(idct2(dtmp) - myidct2(dtmp,2,3)))) 

ttt = zeros(3,3,3);

ttt(1,1,1) = 1; ttt(1,2,1) = 2; ttt(1,3,1) = 3;
ttt(2,1,1) = 1; ttt(2,2,1) = 2; ttt(2,3,1) = 3;
ttt(3,1,1) = 1; ttt(3,2,1) = 2; ttt(3,3,1) = 3;
ttt(1,1,2) = 1; ttt(1,2,2) = 2; ttt(1,3,2) = 3;
ttt(2,1,2) = 1; ttt(2,2,2) = 2; ttt(2,3,2) = 3;
ttt(3,1,2) = 1; ttt(3,2,2) = 2; ttt(3,3,2) = 3;
ttt(1,1,3) = 1; ttt(1,2,3) = 2; ttt(1,3,3) = 3;
ttt(2,1,3) = 1; ttt(2,2,3) = 2; ttt(2,3,3) = 3;
ttt(3,1,3) = 1; ttt(3,2,3) = 2; ttt(3,3,3) = 3;
  
mirt_dctn(ttt)
