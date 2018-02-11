clc

x = (0:100)/100;
T0 = x;
v = z(:,:,1)./z(:,:,2);

for i = Q:-1:2
    k = dsearchn(x',T0')';
    T0 = T0 + v(i,k)/Q;
    sum(x - T0)
end

figure; 
plot(x,T0);