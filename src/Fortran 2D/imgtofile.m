close all 
clear all 

I0 = double(rgb2gray(imread('panda.png')));
I1 = double(rgb2gray(imread('pingouin.png')));

dlmwrite('input/f0.dat',I0)
dlmwrite('input/f1.dat',I1)


F20 = dlmread('results/Transport/f_00020.dat',' ',2); 

F20 = reshape(F20,[],3);

