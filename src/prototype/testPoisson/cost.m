function R = cost(m,f)
    globals;     
    R = 0.5*m.^2./max(f,epsilon);
    
%     R = 0.5*m.^2./f;
%     idx = find(f == 0 & m == 0);
%     R(idx) = 0;
%     idx = find(f < 0);
%     R(idx) = 10^5;

    R = sum(R(:));
end