% Solve -Delta(u) = f
% on ]0,1[
% with u'(0) = u(1) = 0

function u = poisson_neumann(f,N)
    h = 1/N;
    u = zeros(1,N); u1 = u; R = 1.;
    while R > 1e-5
       for i = 2:N-1
           u1(i) = 0.5*(h^2*f(i) + u(i+1) + u(i-1));
       end
       u1(1) = u(2);
       u1(N) = 0;
       R = norm(u-u1)/norm(u)
       u = u1;
    end
end

