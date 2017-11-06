function R = cost(m,f)
% valide seulement en dimension 1 
%ne gere pas le cas infini

R = m.^2./f;
idx = find(f == 0 & m == 0);
R(idx) = 0;

R = sum(sum(R));
end