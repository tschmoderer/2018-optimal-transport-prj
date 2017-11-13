function U = interpolation_adj(V,N,Q)
    m = reshape(V(1:(N+1)*(Q+1)),Q+1,N+1);
    f = reshape(V((N+1)*(Q+1)+1:end),Q+1,N+1);

    %% Matrice d'interpolation %%
    Interpm = [diag(ones(1,N+1));zeros(1,N+1)] + [zeros(1,N+1);diag(ones(1,N+1))];
    Interpm = Interpm/2; % m = mbar*Im
    Interpf = [diag(ones(1,Q+1)) zeros(Q+1,1)] + [zeros(Q+1,1) diag(ones(1,Q+1))];
    Interpf = Interpf/2;

    Interpm_adj = Interpm';
    Interpf_adj = Interpf';
    
    mbar = m*Interpm_adj;
    fbar = Interpf_adj*f;
    
    U = [reshape(mbar,(N+2)*(Q+1),1);reshape(fbar,(N+1)*(Q+2),1)];
end