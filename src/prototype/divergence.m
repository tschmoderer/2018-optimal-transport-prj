function d = divergence(mbar,fbar)
globals;
d = (mbar(:,2:end) - mbar(:,1:end-1)) + (fbar(2:end,:) - fbar(1:end-1,:));
%d = N*(mbar(:,2:end) - mbar(:,1:end-1)) + Q*(fbar(2:end,:) - fbar(1:end-1,:));
end
