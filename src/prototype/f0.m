function f = f0(x)
mu = 0.5;
s = 0.01;
minimal = 0.0001;
f = minimal + exp(-0.5*((x-mu)/s).^2);
f = f/sum(f);
end