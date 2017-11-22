% construit la matrice de l'opérateur frontiere 
% fichier temporaire pour libérer l'espace de travail
% le vecteur U est donné par [mbar1 mbar2 fbar]


globals;
% Opérateur pour f
Bf = blkdiag([diag(ones(1,(N+1)*(Q+1))) zeros((N+1)*(Q+1),(N+1)*(P+1)*Q)],diag(ones(1,(N+1)*(Q+1))));

% Opérateur pour la prmière composante de m
BMmidlle = zeros(2*N,P*(N+2));
for i = 1:2*N
   if mod(i,2) == 1
       BMmidlle(i,1) = 1;
   else 
       BMmidlle(i,end) = 1;
   end
end

Bm1 = [];
for i = 1:Q+1
    Bm1 = blkdiag(Bm1,blkdiag(eye(P+2) , BMmidlle , eye(P+2)));
end

Bm2 = Bm1;
Bm = blkdiag(Bm1,Bm2); % l'opérateur pour m

B = blkdiag(Bm,Bf);


