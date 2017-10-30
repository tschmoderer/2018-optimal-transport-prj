function Iprime=expansion(I)

amin=min(min(I));
amax=max(max(I));
dim=size(I);

A=[amax 1;amin 1];
vec=inv(A)*[1;0];

for(i=1:dim(1))
    for(j=1:dim(2))
        
        
        Iprime(i,j)=vec(1)*I(i,j)+vec(2);
    end
end
