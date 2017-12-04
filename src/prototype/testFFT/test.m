clc
clear all
close all 

%% test d'adjoint de la div
N = 2; P = 1;

mbar = rand(P+1,N+2); fbar = rand(P+2,N+1);
w = rand(P+1,N+1);
wmbar = zeros(P+1,N+2); wfbar = zeros(P+2,N+1);

d = N*(mbar(:,2:end) - mbar(:,1:end-1)) + P*(fbar(2:end,:) - fbar(1:end-1,:));

wmbar(:,1) = -w(:,1); wmbar(:,end) = w(:,end);
wmbar(:,2:end-1) = w(:,1:end-1) - w(:,2:end);

wfbar(1,:) = -w(1,:); wfbar(end,:) = w(end,:); 
wfbar(2:end-1,:) = w(1:end-1,:) - w(2:end,:);

%wmbar = N*wmbar; wfbar = P*wfbar; %% /!\ Provisoire

errDiv = (sum(sum(d.*w)) - sum(sum(mbar.*wmbar)) - sum(sum(fbar.*wfbar)))/sum(sum(d.*w))

%% Interprétation de div en convolution

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

%% Interprétation de div* en convolution
% 2D 
DAdjmbar = -Dm;
mbar1 = [mbar zeros(size(mbar,1),size(DAdjmbar,2)-1);zeros(size(DAdjmbar,1)-1,size(mbar,2)+size(DAdjmbar,2)-1)];
DAdjmbar1 = zeros(size(mbar1));    DAdjmbar1(1:size(DAdjmbar,1),1:size(DAdjmbar,2)) = DAdjmbar;
cM = conv2(w,DAdjmbar,'full');

% 2D 
DAdjfbar = -Df;
fbar1 = [fbar zeros(size(fbar,1),size(DAdjfbar,2)-1);zeros(size(DAdjfbar,1)-1,size(fbar,2)+size(DAdjfbar,2)-1)];
DAdjfbar1 = zeros(size(fbar1)); DAdjfbar1(1:size(DAdjfbar,1),1:size(DAdjfbar,2)) = DAdjfbar;
cF = conv2(w,DAdjfbar,'full');

ecM = sum(sum(abs(wmbar - cM)))
ecF = sum(sum(abs(wfbar - cF)))

%% Opérateur div.div*
conv(Dm,DAdjmbar) % [-1 2 -1]
conv(Df,DAdjfbar) % [-1 ; 2 ; -1]

%% test obtention du masque 2d directement : 
N = 2; P = 2;
w = rand(P+1,N+1); w = [1 2 3;4 5 6;7 8 9]

A = (-w(1:end-2,:) + 2*w(2:end-1,:) - w(3:end,:));
B = (-w(:,1:end-2) + 2*w(:,2:end-1) - w(:,3:end));
ref = A(:,2:end-1) + B(2:end-1,:); 

mW = zeros(3,3); % masque de w.
mW(2,2) = 2; mW(2,[1,3]) = -1; mW([1,3],2) = -1;

w1 = [fbar zeros(size(w,1),size(mW,2)-1);zeros(size(mW,1)-1,size(fbar,2)+size(mW,2)-1)];
mW1 = zeros(size(w1)); mW1(1:size(mW,1),1:size(mW,2)) = mW;

conv2(w,mW,'full')
conv2(w,mW,'same')

fftmW = fft2(mW);
y = zeros(N+1,P+1);
ffty = fft2(y);

ifft2(ffty./fftmW)




