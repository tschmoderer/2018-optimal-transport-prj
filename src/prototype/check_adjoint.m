% check Boundary adjoint 

mX = rand(Q+1,N+2);
fX = rand(Q+2,N+1);

mY = rand(Q+1,2);
fY = rand(2,N+1);

b = boundary(mX,fX);

s1 = sum([sum(b.m.*mY) sum((b.f.*fY)')]);

[mA fA] = boundary_adjoint(mY(:,1),mY(:,2),fY(1,:),fY(2,:),N,Q);

s2 = sum([sum(mX.*mA) sum((fX.*fA)')]);

e1 = abs(s1 - s2)/s1

% check Divergence adjoint

mbarX = rand(Q+1,N+2);
fbarX = rand(Q+2,N+1);

dY = rand(Q+1,N+1);

dX = divergence(mbarX,fbarX);
[mbarY fbarY] = divergence_adjoint(dY);
s1 = sum(sum(dX.*dY));
s2 = sum([sum(mbarX.*mbarY) sum((fbarX.*fbarY)')]);

e2 = abs(s1 - s2)/s1

% check Interpolation adjoint

mbarX = rand(Q+1,N+2);
fbarX = rand(Q+2,N+1);
IYm = rand(Q+1,N+1);
IYf = rand(Q+1,N+1);

[IXm IXf] = interpolation(mbarX,fbarX);
[mbarY fbarY] = interpolation_adjoint(IYm,IYf);

s1 = sum([sum(IXm.*IYm) sum((IXf.*IYf)')]);
s2 = sum([sum(mbarX.*mbarY) sum((fbarX.*fbarY)')]);

e3 = abs(s1 - s2)/s1