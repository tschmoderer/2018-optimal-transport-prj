% construit la matrice de divergence 
% fichier temporaire pour libérer l'espace de travail
% le vecteur U est donné par [mbar1 mbar2 fbar]


globals;
%% calcul divergence en X %%

% bloc élémentaire -- on a une demi grille décalée %
blkdX = zeros((N+1)*(P+2),(N+2)*(P+2));

for i = 1:(N+1)*(P+2)
     for j = 1:(N+2)*(P+2)
         if i == j 
            blkdX(i,j) = -1;
         elseif j == i + P + 2
            blkdX (i,j) = 1;
         end
     end
end
blkdX = N*sparse(blkdX);

% matrice de somme de la demi grille %
blkSdivX = zeros(P+1,P+2);
for i = 1:(P+1)
    for j = 1:(P+2)
        if i == j 
            blkSdivX(i,j) = 1;
        elseif j == i+1
            blkSdivX(i,j) = 1;
        end
    end
end
blkSdivX = sparse(blkSdivX);
SdivX = [];
for i = 1:N+1
    SdivX = blkdiag(SdivX,blkSdivX);
end
SdivX = sparse(SdivX);
% la matrice de divX est la composée des deux %
DmX = SdivX*blkdX;
DmX = sparse(DmX);
clear blkdX SdivX blkSdivX
%% calcul divergence en X %%

% bloc élémentaire -- on a une demi grille décalée %
blkdY = zeros(P+1,P+2);
for i = 1:(P+1)
    for j = 1:(P+2)
        if i == j 
            blkdY(i,j) = -1;
        elseif j == i+1
            blkdY(i,j) = 1;
        end
    end
end
blkdY = P*sparse(blkdY);
blkm2 = [];
for i = 1:N+2
    blkm2 = blkdiag(blkm2,blkdY);
end
blkm2 = sparse(blkm2);
% matrice de somme de la demi grille %
SdivY = zeros((N+1)*(P+1),(N+2)*(P+1));
for i = 1:(N+1)*(P+1)
    for j = 1:(N+2)*(P+1)
        if i == j 
            SdivY(i,j) = 1;
        elseif j == i + P + 1
            SdivY(i,j) = 1;
        end
    end
end
SdivY = sparse(SdivY);
% la matrice de divY est la composée des deux %
DmY = SdivY*blkm2;
DmY = sparse(DmY);
clear SdivY blkm2 blkdY
%% divergence pour la composante m %%
DDX = [];
DDY = [];
for i = 1:Q+1
    DDX = blkdiag(DDX,DmX);
    DDY = blkdiag(DDY,DmY);
end
DDX = sparse(DDX); DDY = sparse(DDY);
Dm = [DDX DDY];
Dm = sparse(Dm);
clear DDX DDY

%% matrice de la divergence en f %
Df = zeros((N+1)*(P+1)*(Q+1),(N+1)*(P+1)*(Q+2));
for i = 1:(N+1)*(P+1)*(Q+1)
    for j = 1:(N+1)*(P+1)*(Q+2)
        if i == j 
            Df(i,j) = -1;
        elseif j == i + N*(P+1) + P +1
            Df(i,j) = 1;
        end 
    end
end
Df = Q*sparse(Df);
%% matrice de divergence %%
D = [Dm Df];
D = sparse(D);
clear Dm Df


