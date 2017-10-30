function normecarree_I=norme_grad(I,dim);

gradI_col=gradient_col(I,dim);
%figure;imagesc(gradI_col),colormap(gray);

gradI_lig=gradient_lig(I,dim);
%figure;imagesc(gradI_lig),colormap(gray);
normecarree_I=gradI_col.^2+gradI_lig.^2;
