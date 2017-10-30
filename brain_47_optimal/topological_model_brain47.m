clear all
close all
clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         SEGMENTATION MODEL UNDER TOPOLOGICAL CONSTRAINTS           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AUTHOR: LE GUYADER CAROLE
% LAST MODIFIED: 07-05-2006
%----------------------------------------------------------------------



%%%%%%%%%%%%%%brain%%%%%%%%%%%%
nu=imread('left_sag_047.tif');
[M N]=size(nu);
nu=[nu;zeros(20,N)];


nu=double(nu);
[M N]=size(nu);
nu=double(nu)/255;

%%%%%%%%%%%%%filtrage prealable de l'image%%%%%%%%%%%
%nu=filtrage_gaussien(nu);%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%figure,imagesc(nu),colormap(gray);

% visualize the data nu in Matlab (rescaled)
% imagesc(nu); axis image; axis equal; axis off; colormap(gray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%Initial 3D surface%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


h=1;
n_lig=[1:M];
n_col=[1:N];
x=h*n_lig';
y=h*n_col;
x0=x*ones(1,N);
y0=ones(M,1)*y;




phi0=-sqrt((x0-(70*h)).^2+...
    (y0-(90*h)).^2)+34;

figure;imagesc(nu);colormap(gray);
hold on;
contour(1:N,1:M,phi0,...
    [0 0],'b');
axis off;

Phi=phi0;
%figure;mesh(Phi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Creation of the edge-detector function%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

normecarree_I=norme_grad(nu,[M N]);
g=1./(1+normecarree_I.^1);
g=expansion(g);
%figure;imagesc(g),colormap(gray);
%colorbar

Phi=geodesic_brain_47(Phi,nu,g);
figure;imagesc(nu);colormap(gray);
hold on;
contour(1:N,1:M,Phi,...
    [0 0],'r');
axis off;