function b = boundary(mbar,fbar)
    b.m = [mbar(1,:);mbar(end,:)];
    b.f = [fbar(:,1),fbar(:,end)];
end 