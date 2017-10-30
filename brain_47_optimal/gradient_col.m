function gradI_col=gradient_col(I,dim)

h=[-1,0,1;-2,0,2;-1,0,1];
gradI_col=conv2(I,h,'same');
column1=I(:,2)-I(:,1);
columndim=I(:,dim(2))-I(:,dim(2)-1);

gradI_col(:,1)=column1;
gradI_col(:,dim(2))=columndim;
