function f = f1(x)
mu = 0.9;
s = 0.005;
f = exp(-0.5*((x-mu)/s).^2)/(s*sqrt(2*pi));
end