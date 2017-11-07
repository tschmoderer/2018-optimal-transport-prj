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
globals;
mbar = rand(Q+1,N+2); fbar = rand(Q+2,N+1);
m = rand(Q+1,N+1); f = rand(Q+1,N+1); 

[Imbar Ifbar] = interpolation(mbar,fbar);
[Iadjm Iadjf] = interpolation_adjoint(m,f); 

s1 = sum(sum(Imbar.*m)) + sum(sum(Ifbar.*f));
s2 = sum(sum(Iadjm.*mbar)) + sum(sum(Iadjf.*fbar));

e3 = abs(s1 - s2)/s1