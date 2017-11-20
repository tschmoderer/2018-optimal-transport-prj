clc
clear all
close


N = 1;
P = 1;
Q = 1;

%% test divergence %%

mbar = rand(P+2,N+2,Q+1,2);
fbar = rand(P+1,N+1,Q+2);

mbar(:,:,1,1) = [3 6 8;3 5 6;7 7 6];
mbar(:,:,2,1) = [6 8 10;6 10 4;1 4 9];
mbar(:,:,1,2) = [3 6 8;3 5 6;7 7 6];
mbar(:,:,2,2) = [6 8 10;6 10 4;1 4 9];

% reference %
dX = mbar(:,2:end,:,1) - mbar(:,1:end-1,:,1); % div selon x
dY = mbar(2:end,:,:,2) - mbar(1:end-1,:,:,2); % div selon y
d = dX(2:end,:,:) + dX(1:end-1,:,:) + dY(:,2:end,:) + dY(:,1:end-1,:) + fbar(:,:,2:end) - fbar(:,:,1:end-1); % div totale

U = reshape([mbar(:,:,:,1) ; mbar(:,:,:,2)],[],1);


% divergence en X on agit mbar1
blkm1 = zeros((N+1)*(P+2),(N+2)*(P+2));
 
for i = 1:(N+1)*(P+2)
     for j = 1:(N+2)*(P+2)
         if i == j 
            blkm1(i,j) = -1;
         elseif j == i + P + 2
            blkm1 (i,j) = 1;
         end
     end
end

blkSdivX = zeros(P+1,P+2);
for i = 1:(P+1)
    for j = 1:(P+2)
        if i == j 
            blkSdivX(i,j) = 1;
        elseif j == i+1
            blkSdivX(i,j) = 1;
        end
    end
end
SdivX = [];
for i = 1:N+1
    SdivX = blkdiag(SdivX,blkSdivX);
end

DmX = SdivX*blkm1;

% divergence en Y on agit sur mbar2

blkdmY = zeros(P+1,P+2);
for i = 1:(P+1)
    for j = 1:(P+2)
        if i == j 
            blkdmY(i,j) = -1;
        elseif j == i+1
            blkdmY(i,j) = 1;
        end
    end
end
blkm2 = [];
for i = 1:N+2
    blkm2 = blkdiag(blkm2,blkdmY);
end

SdivY = zeros((N+1)*(P+1),(N+2)*(P+1));
for i = 1:(N+1)*(P+1)
    for j = 1:(N+2)*(P+1)
        if i == j 
            SdivY(i,j) = 1;
        elseif j == i + P + 1
            SdivY(i,j) = 1;
        end
    end
end

DmY = SdivY*blkm2;

Dm = [];
for i = 1:Q+1
   Dm = blkdiag(Dm,[DmX DmY]); 
end

% a reprendre
Df = zeros((N+1)*(P+1)*(Q+1),(N+1)*(P+1)*(Q+2));
for i = 1:(N+1)*(P+1)*(Q+1)
    for j = 1:(N+1)*(P+1)*(Q+2)
        if i ==j 
            Df(i,j) = -1;
        elseif j == i + P + 1
            Df(i,j) = -1;
        end 
    end
end

% dX(2:end,:,:) + dX(1:end-1,:,:) + dY(:,2:end,:) + dY(:,1:end-1,:)
% reshape([DmX DmY]*U,P+1,N+1,Q+1)














% DmX = [];
% for i = 1:Q+1
%     DmX = blkdiag(DmX,blkm1);
% end
% 
% sum(sum(sum(abs(reshape(DmX*reshape(mbar(:,:,:,1),[],1),[P+1,N+2,Q+1]) - dX))))
% 
% 
% 
% blkY = zeros(P+1,P+2);
% for i = 1:P+1
%     for j = 1:P+2
%         if i == j 
%             blkY(i,j) = -1;
%         elseif j == i+1
%             blkY(i,j) =1;
%         end
%     end
% end
%  
% DmY = [];
% for i = 1:(N+2)*(Q+1)
%     DmY = blkdiag(DmY,blkY);
% end
% 
% sum(sum(sum(abs(reshape(DmY*reshape(mbar(:,:,:,2),[],1),[N+1,P+2,Q+1]) - dY))))

%Df*reshape(fbar,[],1);
%D = [N*Dm1 P*Dm2 Q*Df];


%ù fin test div %%