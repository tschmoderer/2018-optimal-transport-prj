function [m f] = interpolation(mbar,fbar)

m = mbar(:,2:end) + mbar(:,1:end-1);
f = fbar(2:end,:) + fbar(1:end-1,:);

m = 0.5*m;
f = 0.5*f;

%I = zeros(7,6) + [diag(ones(1,6));zeros(1,6)] + [zeros(1,6);diag(ones(1,6))]
end