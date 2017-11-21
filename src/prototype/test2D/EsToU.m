function U = EsToU(mbar,fbar)
    U = [reshape(mbar(:,:,:,1),[],1) ; reshape(mbar(:,:,:,2),[],1) ; reshape(fbar,[],1)];
end

