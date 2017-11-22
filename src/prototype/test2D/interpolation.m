% construit la matrice d'interpolation 
% fichier temporaire pour libérer l'espace de travail
% le vecteur U est donné par [mbar1 mbar2 fbar]

globals;
%% matrice de l'interpolation en f %%
If = zeros((N+1)*(P+1)*(Q+1),(N+1)*(P+1)*(Q+2));
for i = 1:(N+1)*(P+1)*(Q+1)
    for j = 1:(N+1)*(P+1)*(Q+2)
        if i == j 
            If(i,j) = 1;
        elseif j == i + N*(P+1) + P +1
            If(i,j) = 1;
        end 
    end
end

%% matrice d'interpolation en m1 %%

blk1Im1 = zeros((N+1)*(P+2),(N+2)*(P+2));

for i = 1:(N+1)*(P+2)
    for j = 1:(N+2)*(P+2)
        if i == j
            blk1Im1(i,j) = 1;
        elseif j == i + P + 2
            blk1Im1(i,j) = 1;
        end
    end
end


blk = zeros(P+1,P+2); 
for i = 1:P+1
    for j = 1:P+2
        if i == j 
            blk(i,j) = 1;
        elseif j == i+1
            blk(i,j) = 1;
        end
    end
end

blk2Im1 = [];

for i = 1:(N+1)
    blk2Im1 = blkdiag(blk2Im1,blk);
end

Im1 = [];

for i = 1:Q+1;
    Im1 = blkdiag(Im1,blk2Im1*blk1Im1);
end

I = blkdiag(0.25*Im1,0.25*Im1,0.5*If);
