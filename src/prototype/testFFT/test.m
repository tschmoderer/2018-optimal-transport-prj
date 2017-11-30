clc
clear all
close all 

N = 2;
Q = 1;

m = rand(Q+1,N+1);
f = rand(Q+1,N+1); 

mbar = rand(Q+1,N+2);
fbar = rand(Q+2,N+1);

%% Test pour faire de la divergence un opérateur de convolution %%

% reference %
dm = mbar(:,2:end) - mbar(:,1:end-1);
df = fbar(2:end,:) - fbar(1:end-1,:);

% 2D 
Dm = [1 -1];
mbar1 = [mbar zeros(size(mbar,1),size(Dm,2)-1);zeros(size(Dm,1)-1,size(mbar,2)+size(Dm,2)-1)];
Dm1 = zeros(size(mbar1));    Dm1(1:size(Dm,1),1:size(Dm,2)) = Dm;
c1 = ifft2(fft2(mbar1).*fft2(Dm1));
tm = conv2(mbar,Dm,'full');

% 2D 
Df = [1 ; -1];
fbar1 = [fbar zeros(size(fbar,1),size(Df,2)-1);zeros(size(Df,1)-1,size(fbar,2)+size(Df,2)-1)];
Df1 = zeros(size(fbar1)); Df1(1:size(Df,1),1:size(Df,2)) = Df;
c1 = ifft2(fft2(fbar1).*fft2(Df1));
tf = conv2(fbar,Df,'full'); 

eM = sum(sum(abs(dm-tm(:,2:end-1))))
eF = sum(sum(abs(df-tf(2:end-1,:))))

% essaie de faire Dm***mbar = dm // *** = conv
Dm = [1 -1];
mbar1 = [mbar(:,2:end-1) zeros(size(mbar(:,2:end-1),1),size(Dm,2)-1);zeros(size(Dm,1)-1,size(mbar(:,2:end-1),2)+size(Dm,2)-1)];
Dm1 = zeros(size(mbar1));    Dm1(1:size(Dm,1),1:size(Dm,2)) = Dm;

tm2 = conv2(mbar(:,2:end-1),Dm)
%dm
%tm

