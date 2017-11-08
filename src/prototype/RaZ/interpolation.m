function [m,f] = interpolation(mbar,fbar)
    m = (mbar(1:end-1,:) + mbar(2:end,:))/2;
    f = (fbar(:,1:end-1) + fbar(:,2:end))/2;
end