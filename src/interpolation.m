function V = interpolation(U)
    globals;
%     mbar = reshape(U(1:(Q+1)*(N+2)),Q+1,N+2);
%     fbar = reshape(U((Q+1)*(N+2)+1:end),Q+2,N+1);
%     
%     %% Matrice d'interpolation %%
%     Interpm = [diag(ones(1,N+1));zeros(1,N+1)] + [zeros(1,N+1);diag(ones(1,N+1))];
%     Interpm = Interpm/2; % m = mbar*Im
%     Interpf = [diag(ones(1,Q+1)) zeros(Q+1,1)] + [zeros(Q+1,1) diag(ones(1,Q+1))];
%     Interpf = Interpf/2;
%     
%     m = mbar*Interpm;
%     f = Interpf*fbar;
%     
%     V = [reshape(m,(N+1)*(Q+1),1);reshape(f,(N+1)*(Q+1),1)];
    V = Interp*U;
end