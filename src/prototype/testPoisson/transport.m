clc
clear all
close all

%% Initialisation %%

N = 29;
Q = 31;


alpha = 1.998; g = 0.0125;

[XX,YY] = meshgrid(linspace(0,1,N+1),linspace(0,1,Q+1)); YY = flipud(YY);

wmbar0 = zeros(Q+1,N+2); wmbar1 = zeros(Q+1,N+2); zmbar = zeros(Q+1,N+2);
wfbar0 = zeros(Q+2,N+1); wfbar1 = zeros(Q+2,N+1); zfbar = zeros(Q+2,N+1);
wm0 = zeros(Q+1,N+1); wm1 = zeros(Q+1,N+1); zm = zeros(Q+1,N+1);
wf0 = zeros(Q+1,N+1); wf1 = zeros(Q+1,N+1); zf = zeros(Q+1,N+1);


% Itérations
niter = 1000;
cout = zeros(1,niter);
div = zeros(1,niter);
minF = zeros(1,niter);

tic;
for l = 1:niter
    % etape 1
    [wmbar1,wfbar1,wm1,wf] = proxG1(2*zmbar-wmbar0,2*zfbar-wfbar0,2*zm-wm0,2*zf-wf0,N,Q);    
    
    wmbar1 = wmbar0 + alpha*(wmbar1-zmbar); wfbar1 = wfbar0 + alpha*(wfbar1-zfbar);
    wm1 = wm0 + alpha*(wm1-zm); wf1 = wf0 + alpha*(wf1 - zf);
    
    % etape 2
    [zmbar,zfbar,zm,zf] = proxG2(wmbar1,wfbar1,wm1,wf1,N,Q);
   
    % iterate
    wmbar0 = wmbar1; wfbar0 = wfbar1;
    wm0 = wm1; wf0 = wf1;
    
    
    surf(XX,YY,zf)
    xlabel('x');
    ylabel('t');
    zlabel('f');
    title(['itération : ',num2str(l)]);
    drawnow

    
    cout(l) = cost(zm,zf);
   % div(l)  = sum(D*zU0);
    minF(l) = min(zf(:));
end
toc

figure;
subplot(3,1,1)
plot([1:niter],cout);
title('cout')
subplot(3,1,2)
plot([1:niter],div)
title('div')
subplot(3,1,3)
plot([1:niter],minF);
title('min de f');
