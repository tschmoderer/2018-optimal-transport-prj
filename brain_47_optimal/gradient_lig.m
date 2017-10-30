function gradI_lig=gradient_lig(I,dim)

h=[-1,-2,-1;0,0,0;1,2,1];
gradI_lig=conv2(I,h,'same');

ligne1=I(2,:)-I(1,:);
lignedim=I(dim(1),:)-I(dim(1)-1,:);

gradI_lig(1,:)=ligne1;
gradI_lig(dim(1),:)=lignedim;
