function res = poisson3d_Neumann(f,lx,ly,lz)

%
%
%  res = poisson3d_Neumann(f,lx,ly,lz)
%
%


if(exist('lx','var')==0)
  lx                  = 1;
end
if(exist('ly','var')==0)
  ly                  = 1;
end
if(exist('lz','var')==0)
  lz                  = 1;
end
N                     = size(f,1);
M                     = size(f,2);
R                     = size(f,3);
hx                    = lx/N;
hy                    = ly/M;
hz                    = lz/R;

dn                    = 0:1:N-1;
depn                  = 2*cos(pi*dn/N)-2;
dm                    = 0:1:M-1;
depm                  = 2*cos(pi*dm/M)-2;
dr                    = 0:1:R-1;
depr                  = 2*cos(pi*dr/R)-2;

denom2                = repmat(depn(:)/hx^2,[1,M,R]) + repmat(depm(:)'/hy^2,[N,1,R]) + permute(repmat(depr(:)/hz^2,[1,N,M]),[2 3 1]);
denom2(denom2(:)==0)  = 1;

%for i=1:N
%    for j=1:M
%       for k=1:R
%            denom(i,j,k) = depn(i)/hx^2 + depm(j)/hy^2 + depr(k)/hz^2;
%            if (denom(i,j,k)==0)
%                denom(i,j,k) = 1;
%            end
%       end
%    end
%end
%max(abs(denom(:)-denom2(:)))


fhat                  = mirt_dctn(f);
uhat                  = -(fhat)./denom2;
res                   = mirt_idctn(uhat);

