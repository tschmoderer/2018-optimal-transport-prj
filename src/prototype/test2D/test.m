clc
clear all
close

%%% test reshape 3D %%%

N = 1;
Q = 1;

m = zeros(Q+1,N+1,2);
m(:,:,1) = [1 2 ; 3 4]
m(:,:,2) = [5 6 ; 7 8]