  #include "mex.h"
  #include <stdio.h>
  #include <math.h>
  #include <stdlib.h>


  /***********************************/
  /*********global variables*********/
  /**********************************/
  /*********time discretization*******/
   double dt=0.5;
  /*********space discretization********/
   double  h=1.0;
  /******parameter balancing the constant motion*******/
   double regul=-0.2;
   /********number of iterations in the reinitialization process***/
   int compteur=20;
   /********parameter epsilon in the regularized versions of H and delta********/
   double epsilon=1.0;
   /********parameter balancing the topological constraint*************/
   double mu=-0.2;


   void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]){

   int Itermax=400;
   int i,j,k;
   int mrows,ncols;
   int level=1;
   int window=2*level;
   double energy=0.0;
   double ** nu;
   double ** Phi;
   double ** sgnPhi0;
   double ** Phik;
   double ** Phix;
   double ** Phiy;
   double ** Phiz;
   double ** Psi;
   double ** contrib;
   double ** g;
   double **f;
   double ** topo;
   double ** norme_Phi;
   double * flig_balloon;
   double * fcol_balloon;
   double * topolig;
   double * topocol;
   double * smlig;
   double * smcol;
   double * sol_lig;
   double * sol_col;
   double * kcpy;
   double * v_Phi;
   double * v_Phi_col;
   double * v_Phi_lig;
   double ** A_x;
   double ** A_y;
   double ** mat_col;
   /**************************************************************************/
   /*************declaration of the functions used in the main program********/
   /**************************************************************************/

   void balloon_force(int ,int ,double ** ,double ** ,double ** );
   void norme_gradPhi(int ,int ,double ** ,double ** );
   void vectorisation_colonne(int,int,double**,double *);
   void vectorisation_ligne(int,int,double**,double *);
   void matrix_A_x2(int,int,double ** ,double ** ,double ** );
   void matrix_A_y2(int,int,double ** ,double ** ,double ** );
   void thomas(int ,double ** ,double * ,double *);
   void matrix_col(int ,int ,double * ,double** );
   void matrix_lig(int ,int ,double * ,double** );
   void reinitialization(int , int ,double **,double **,double **);
   void topological_force(int,int,double **,double **,int,int,double ** ,
                             double ** ,double ** , double ** );
   double heaviside(double,double);
   double dirac(double,double);

   if(nrhs!=3){
   mexErrMsgTxt("Two inputs required.");
   }
   if(nlhs>1){
   mexErrMsgTxt("One output argument.");
   }

   mrows=mxGetM(prhs[0]);
   ncols=mxGetN(prhs[0]);

   nu=(double**)malloc(mrows*sizeof(double *));
   Phi=(double**)malloc(mrows*sizeof(double *));
   Phik=(double**)malloc(mrows*sizeof(double *));
   sgnPhi0=(double**)malloc(mrows*sizeof(double *));
   Phix=(double**)malloc((mrows+2*window)*sizeof(double *));
   Phiy=(double**)malloc((mrows+2*window)*sizeof(double *));
   Phiz=(double**)malloc((mrows+2*window)*sizeof(double *));
   Psi=(double**)malloc((mrows+2*window)*sizeof(double *));
   contrib=(double**)malloc((mrows+2*window)*sizeof(double *));
   g=(double**)malloc(mrows*sizeof(double*));
   f=(double**)malloc(mrows*sizeof(double*));
   topo=(double**)malloc(mrows*sizeof(double*));
   norme_Phi=(double**)malloc(mrows*sizeof(double*));
   mat_col=(double**)malloc(mrows*sizeof(double*));
   flig_balloon=(double*)malloc(sizeof(double)*mrows*ncols);
   fcol_balloon=(double*)malloc(sizeof(double)*mrows*ncols);
   topolig=(double*)malloc(sizeof(double)*mrows*ncols);
   topocol=(double*)malloc(sizeof(double)*mrows*ncols);
   v_Phi=(double*)malloc(sizeof(double)*mrows*ncols);
   v_Phi_col=(double*)malloc(sizeof(double)*mrows*ncols);
   v_Phi_lig=(double*)malloc(sizeof(double)*mrows*ncols);
   A_x=(double**)calloc(mrows*ncols,sizeof(double *));
   A_y=(double**)calloc(mrows*ncols,sizeof(double *));
   smlig=(double*)malloc(sizeof(double)*mrows*ncols);
   smcol=(double*)malloc(sizeof(double)*mrows*ncols);
   sol_lig=(double*)malloc(sizeof(double)*mrows*ncols);
   sol_col=(double*)malloc(sizeof(double)*mrows*ncols);

   for(i=0;i<mrows;i++){
   nu[i]=(double*)calloc(ncols,sizeof(double));
   Phi[i]=(double*)calloc(ncols,sizeof(double));
   Phik[i]=(double*)calloc(ncols,sizeof(double));
   sgnPhi0[i]=(double*)calloc(ncols,sizeof(double));
   g[i]=(double*)calloc(ncols,sizeof(double));
   f[i]=(double*)calloc(ncols,sizeof(double));
   norme_Phi[i]=(double*)calloc(ncols,sizeof(double));
   mat_col[i]=(double*)calloc(ncols,sizeof(double));
   Phik[i]=(double*)calloc(ncols,sizeof(double));
   topo[i]=(double*)calloc(ncols,sizeof(double));
   }

   for(i=0;i<mrows+2*window;i++){
   Phix[i]=(double*)calloc(ncols+2*window,sizeof(double));
   Phiy[i]=(double*)calloc(ncols+2*window,sizeof(double));
   Phiz[i]=(double*)calloc(ncols+2*window,sizeof(double));
   Psi[i]=(double*)calloc(ncols+2*window,sizeof(double));
   contrib[i]=(double*)calloc(ncols+2*window,sizeof(double));
   }

   for(i=0;i<mrows*ncols;i++){
   A_x[i]=(double*)calloc(3,sizeof(double));
   A_y[i]=(double*)calloc(3,sizeof(double));
   }

   for(i=0;i<ncols; i++){
            for(j=0; j<mrows; j++){
   Phi[j][i]=(mxGetPr(prhs[0]))[i*mrows+j];
   nu[j][i]=(mxGetPr(prhs[1]))[i*mrows+j];
    g[j][i]=(mxGetPr(prhs[2]))[i*mrows+j];
   }
   }


  norme_gradPhi(mrows,ncols,Phi,norme_Phi);
  balloon_force(mrows,ncols,Phi,g,f);
  topological_force(mrows,ncols,Phi,topo,level,window,Phix,Phiy,Psi,contrib);
  vectorisation_ligne(mrows,ncols,f,flig_balloon);
  vectorisation_colonne(mrows,ncols,f,fcol_balloon);
  vectorisation_ligne(mrows,ncols,Phi,v_Phi_lig);
  vectorisation_colonne(mrows,ncols,Phi,v_Phi_col);
  vectorisation_ligne(mrows,ncols,topo,topolig);
  vectorisation_colonne(mrows,ncols,topo,topocol);

  for(k=0;k<Itermax;k++){

  matrix_A_x2(mrows,ncols,g,norme_Phi,A_x);
  matrix_A_y2(mrows,ncols,g,norme_Phi,A_y);

  for(i=0;i<mrows*ncols;i++){
  smlig[i]=v_Phi_lig[i]+flig_balloon[i]+topolig[i];
  smcol[i]=v_Phi_col[i]+fcol_balloon[i]+topocol[i];
  }

  thomas(mrows*ncols,A_x,smlig,sol_lig);
  thomas(mrows*ncols,A_y,smcol,sol_col);
  matrix_col(mrows,ncols,sol_col,mat_col);
  vectorisation_ligne(mrows,ncols,mat_col,sol_col);

  for(i=0;i<mrows*ncols;i++){
  v_Phi[i]=0.5*(sol_lig[i]+sol_col[i]);

  }
  matrix_lig(mrows,ncols,v_Phi,Phi);

  reinitialization(mrows,ncols,Phi,Phik,sgnPhi0);
  for(i=0;i<mrows;i++){
    for(j=0;j<ncols;j++){
    Phi[i][j]=Phik[i][j];
    }
  }


  norme_gradPhi(mrows,ncols,Phi,norme_Phi);
  balloon_force(mrows,ncols,Phi,g,f);
  topological_force(mrows,ncols,Phi,topo,level,window,Phix,Phiy,Psi,contrib);
  vectorisation_ligne(mrows,ncols,f,flig_balloon);
  vectorisation_colonne(mrows,ncols,f,fcol_balloon);
  vectorisation_ligne(mrows,ncols,Phi,v_Phi_lig);
  vectorisation_colonne(mrows,ncols,Phi,v_Phi_col);
  vectorisation_ligne(mrows,ncols,topo,topolig);
  vectorisation_colonne(mrows,ncols,topo,topocol);


  
  }


  plhs[0]=mxCreateDoubleMatrix(mrows,ncols,mxREAL);
  kcpy=mxGetPr(plhs[0]);

  for(i=0;i<ncols;i++){
            for(j=0;j<mrows;j++){

              kcpy[i*mrows+j]=Phi[j][i];
  }
  }

  for(i=0;i<mrows;i++){
  free(nu[i]);
  free(Phi[i]);
  free(g[i]);
  free(f[i]);
  free(norme_Phi[i]);
  free(mat_col[i]);
  free(Phik[i]);
  free(sgnPhi0[i]);
  }
  for(i=0;i<mrows*ncols;i++){
  free(A_x[i]);
  free(A_y[i]);
  }
  for(i=0;i<mrows+2*window;i++){
  free(Phix[i]);
  free(Phiy[i]);
  free(Phiz[i]);
  free(Psi[i]);
  free(contrib[i]);
  }

  free(nu),free(Phi),free(g),free(f),free(norme_Phi),free(mat_col),free(Phik),free(sgnPhi0),free(A_x),free(A_y);
  free(flig_balloon),free(fcol_balloon),free(v_Phi),free(v_Phi_col),free(v_Phi_lig),free(smlig),free(smcol),free(sol_lig),free(sol_col);
  free(topolig),free(topocol),free(Phix),free(Phiy),free(Phiz),free(contrib),free(Psi);
}







   /**********************************************************/
   /**************definition of the max function**************/
   /**********************************************************/

   double getmax(double a,double b)

  {

   if(a<=b){
     return(b);
   }
   else{return(a);
  }
  }

   /**********************************************************/
   /**************definition of the min function**************/
   /**********************************************************/
  double getmin(double a, double b)
  
  {
  if(a<=b){
     return(a);
   }
   else{return(b);
  }
  }
  
  
   /*************************************************************/
   /************definition of the abs function*******************/
   /*************************************************************/
   
   double getabs(double a){
   if(a>=0){
     return(a);
   }
   else return(-a);
   }
   /*************************************************************/
   /***********definition of the function balloon force**********/
   /*discretization of the gradient of Phi satisfying the entropy*/
   /************************* condition.**************************/
   /*************************************************************/

   void balloon_force(int nblig,int nbcol,double ** Phi,double ** g,double ** f) {

   int i,j;
   double ** norme_Phi;

   norme_Phi=(double**)malloc(nblig*sizeof(double *));
   for(i=0;i<nblig;i++){
   norme_Phi[i]=(double*)calloc(nbcol,sizeof(double));
   }
   /**************************************************/
   /*************upper left corner********************/
   /**************************************************/

   if (regul <=0.0){
   norme_Phi[0][0]=sqrt(getmax(Phi[0][0]-Phi[0][1],0.0)*getmax(Phi[0][0]-Phi[0][1],0.0)+
                        getmin(Phi[0][1]-Phi[0][0],0.0)*getmin(Phi[0][1]-Phi[0][0],0.0)+
                        getmax(Phi[0][0]-Phi[1][0],0.0)*getmax(Phi[0][0]-Phi[1][0],0.0)+
                        getmin(Phi[1][0]-Phi[0][0],0.0)*getmin(Phi[1][0]-Phi[0][0],0.0));

   }
   else{
   norme_Phi[0][0]=sqrt(getmin(Phi[0][0]-Phi[0][1],0.0)*getmin(Phi[0][0]-Phi[0][1],0.0)+
                        getmax(Phi[0][1]-Phi[0][0],0.0)*getmax(Phi[0][1]-Phi[0][0],0.0)+
                        getmin(Phi[0][0]-Phi[1][0],0.0)*getmin(Phi[0][0]-Phi[1][0],0.0)+
                        getmax(Phi[1][0]-Phi[0][0],0.0)*getmax(Phi[1][0]-Phi[0][0],0.0));
   }
   
   /**************************************************/
   /*************upper right corner*******************/
   /**************************************************/
   
   if (regul <=0.0){
   norme_Phi[0][nbcol-1]=sqrt(getmax(Phi[0][nbcol-1]-Phi[0][nbcol-2],0.0)*getmax(Phi[0][nbcol-1]-Phi[0][nbcol-2],0.0)+
                        getmin(Phi[0][nbcol-2]-Phi[0][nbcol-1],0.0)*getmin(Phi[0][nbcol-2]-Phi[0][nbcol-1],0.0)+
                        getmax(Phi[0][nbcol-1]-Phi[1][nbcol-1],0.0)*getmax(Phi[0][nbcol-1]-Phi[1][nbcol-1],0.0)+
                        getmin(Phi[1][nbcol-1]-Phi[0][nbcol-1],0.0)*getmin(Phi[1][nbcol-1]-Phi[0][nbcol-1],0.0));

   }
   else{
   norme_Phi[0][nbcol-1]=sqrt(getmin(Phi[0][nbcol-1]-Phi[0][nbcol-2],0.0)*getmin(Phi[0][nbcol-1]-Phi[0][nbcol-2],0.0)+
                        getmax(Phi[0][nbcol-2]-Phi[0][nbcol-1],0.0)*getmax(Phi[0][nbcol-2]-Phi[0][nbcol-1],0.0)+
                        getmin(Phi[0][nbcol-1]-Phi[1][nbcol-1],0.0)*getmin(Phi[0][nbcol-1]-Phi[1][nbcol-1],0.0)+
                        getmax(Phi[1][nbcol-1]-Phi[0][nbcol-1],0.0)*getmax(Phi[1][nbcol-1]-Phi[0][nbcol-1],0.0));
   }

   /*****************************************************/
   /*****************lower left corner*******************/
   /*****************************************************/
   
   if (regul <=0.0){
   norme_Phi[nblig-1][0]=sqrt(getmax(Phi[nblig-1][0]-Phi[nblig-1][1],0.0)*getmax(Phi[nblig-1][0]-Phi[nblig-1][1],0.0)+
                        getmin(Phi[nblig-1][1]-Phi[nblig-1][0],0.0)*getmin(Phi[nblig-1][1]-Phi[nblig-1][0],0.0)+
                        getmax(Phi[nblig-1][0]-Phi[nblig-2][0],0.0)*getmax(Phi[nblig-1][0]-Phi[nblig-2][0],0.0)+
                        getmin(Phi[nblig-2][0]-Phi[nblig-1][0],0.0)*getmin(Phi[nblig-2][0]-Phi[nblig-1][0],0.0));

   }

   else{
   
   norme_Phi[nblig-1][0]=sqrt(getmin(Phi[nblig-1][0]-Phi[nblig-1][1],0.0)*getmin(Phi[nblig-1][0]-Phi[nblig-1][1],0.0)+
                        getmax(Phi[nblig-1][1]-Phi[nblig-1][0],0.0)*getmax(Phi[nblig-1][1]-Phi[nblig-1][0],0.0)+
                        getmin(Phi[nblig-1][0]-Phi[nblig-2][0],0.0)*getmin(Phi[nblig-1][0]-Phi[nblig-2][0],0.0)+
                        getmax(Phi[nblig-2][0]-Phi[nblig-1][0],0.0)*getmax(Phi[nblig-2][0]-Phi[nblig-1][0],0.0));


   }


   /*****************************************************/
   /*****************lower right corner******************/
   /*****************************************************/
   
   if (regul <=0.0){
   norme_Phi[nblig-1][nbcol-1]=sqrt(getmax(Phi[nblig-1][nbcol-1]-Phi[nblig-1][nbcol-2],0.0)*getmax(Phi[nblig-1][nbcol-1]-Phi[nblig-1][nbcol-2],0.0)+
                        getmin(Phi[nblig-1][nbcol-2]-Phi[nblig-1][nbcol-1],0.0)*getmin(Phi[nblig-1][nbcol-2]-Phi[nblig-1][nbcol-1],0.0)+
                        getmax(Phi[nblig-1][nbcol-1]-Phi[nblig-2][nbcol-1],0.0)*getmax(Phi[nblig-1][nbcol-1]-Phi[nblig-2][nbcol-1],0.0)+
                        getmin(Phi[nblig-2][nbcol-1]-Phi[nblig-1][nbcol-1],0.0)*getmin(Phi[nblig-2][nbcol-1]-Phi[nblig-1][nbcol-1],0.0));

   }

   else{
   norme_Phi[nblig-1][nbcol-1]=sqrt(getmin(Phi[nblig-1][nbcol-1]-Phi[nblig-1][nbcol-2],0.0)*getmin(Phi[nblig-1][nbcol-1]-Phi[nblig-1][nbcol-2],0.0)+
                        getmax(Phi[nblig-1][nbcol-2]-Phi[nblig-1][nbcol-1],0.0)*getmax(Phi[nblig-1][nbcol-2]-Phi[nblig-1][nbcol-1],0.0)+
                        getmin(Phi[nblig-1][nbcol-1]-Phi[nblig-2][nbcol-1],0.0)*getmin(Phi[nblig-1][nbcol-1]-Phi[nblig-2][nbcol-1],0.0)+
                        getmax(Phi[nblig-2][nbcol-1]-Phi[nblig-1][nbcol-1],0.0)*getmax(Phi[nblig-2][nbcol-1]-Phi[nblig-1][nbcol-1],0.0));

   
   }
   /******************************************************/
   /***************left edge without the corners**********/
   /******************************************************/
   
   if (regul <=0.0){

   for(i=1;i<nblig-1;i++){

   norme_Phi[i][0]=sqrt(getmax(Phi[i][0]-Phi[i][1],0.0)*getmax(Phi[i][0]-Phi[i][1],0.0)+
                        getmin(Phi[i][1]-Phi[i][0],0.0)*getmin(Phi[i][1]-Phi[i][0],0.0)+
                        getmax(Phi[i][0]-Phi[i-1][0],0.0)*getmax(Phi[i][0]-Phi[i-1][0],0.0)+
                        getmin(Phi[i+1][0]-Phi[i][0],0.0)*getmin(Phi[i+1][0]-Phi[i][0],0.0));

   }
   }
   else{
   
   for(i=1;i<nblig-1;i++){

   norme_Phi[i][0]=sqrt(getmin(Phi[i][0]-Phi[i][1],0.0)*getmin(Phi[i][0]-Phi[i][1],0.0)+
                        getmax(Phi[i][1]-Phi[i][0],0.0)*getmax(Phi[i][1]-Phi[i][0],0.0)+
                        getmin(Phi[i][0]-Phi[i-1][0],0.0)*getmin(Phi[i][0]-Phi[i-1][0],0.0)+
                        getmax(Phi[i+1][0]-Phi[i][0],0.0)*getmax(Phi[i+1][0]-Phi[i][0],0.0));
   }
   }


   /***********************************************************/
   /****************right edge without the corners*************/
   /***********************************************************/

   if (regul <=0.0){

   for(i=1;i<nblig-1;i++){

   norme_Phi[i][nbcol-1]=sqrt(getmax(Phi[i][nbcol-1]-Phi[i][nbcol-2],0.0)*getmax(Phi[i][nbcol-1]-Phi[i][nbcol-2],0.0)+
                        getmin(Phi[i][nbcol-2]-Phi[i][nbcol-1],0.0)*getmin(Phi[i][nbcol-2]-Phi[i][nbcol-1],0.0)+
                        getmax(Phi[i][nbcol-1]-Phi[i-1][nbcol-1],0.0)*getmax(Phi[i][nbcol-1]-Phi[i-1][nbcol-1],0.0)+
                        getmin(Phi[i+1][nbcol-1]-Phi[i][nbcol-1],0.0)*getmin(Phi[i+1][nbcol-1]-Phi[i][nbcol-1],0.0));

   }
   }

   else{
   for(i=1;i<nblig-1;i++){
   norme_Phi[i][nbcol-1]=sqrt(getmin(Phi[i][nbcol-1]-Phi[i][nbcol-2],0.0)*getmin(Phi[i][nbcol-1]-Phi[i][nbcol-2],0.0)+
                        getmax(Phi[i][nbcol-2]-Phi[i][nbcol-1],0.0)*getmax(Phi[i][nbcol-2]-Phi[i][nbcol-1],0.0)+
                        getmin(Phi[i][nbcol-1]-Phi[i-1][nbcol-1],0.0)*getmin(Phi[i][nbcol-1]-Phi[i-1][nbcol-1],0.0)+
                        getmax(Phi[i+1][nbcol-1]-Phi[i][nbcol-1],0.0)*getmax(Phi[i+1][nbcol-1]-Phi[i][nbcol-1],0.0));

   }
   }

   /***********************************************************/
   /****************upper edge without the corners*************/
   /***********************************************************/

   if (regul <=0.0){

   for(j=1;j<nbcol-1;j++){

   norme_Phi[0][j]=sqrt(getmax(Phi[0][j]-Phi[0][j-1],0.0)*getmax(Phi[0][j]-Phi[0][j-1],0.0)+
                        getmin(Phi[0][j+1]-Phi[0][j],0.0)*getmin(Phi[0][j+1]-Phi[0][j],0.0)+
                        getmax(Phi[0][j]-Phi[1][j],0.0)*getmax(Phi[0][j]-Phi[1][j],0.0)+
                        getmin(Phi[1][j]-Phi[0][j],0.0)*getmin(Phi[1][j]-Phi[0][j],0.0));

   }
   }
   
   else{
   for(j=1;j<nbcol-1;j++){

   norme_Phi[0][j]=sqrt(getmin(Phi[0][j]-Phi[0][j-1],0.0)*getmin(Phi[0][j]-Phi[0][j-1],0.0)+
                        getmax(Phi[0][j+1]-Phi[0][j],0.0)*getmax(Phi[0][j+1]-Phi[0][j],0.0)+
                        getmin(Phi[0][j]-Phi[1][j],0.0)*getmin(Phi[0][j]-Phi[1][j],0.0)+
                        getmax(Phi[1][j]-Phi[0][j],0.0)*getmax(Phi[1][j]-Phi[0][j],0.0));

   }

   }
   
   /************************************************************/
   /*****************lower edge without the corners*************/
   /************************************************************/
   
   if (regul <=0.0){

   for(j=1;j<nbcol-1;j++){

   norme_Phi[nblig-1][j]=sqrt(getmax(Phi[nblig-1][j]-Phi[nblig-1][j-1],0.0)*getmax(Phi[nblig-1][j]-Phi[nblig-1][j-1],0.0)+
                        getmin(Phi[nblig-1][j+1]-Phi[nblig-1][j],0.0)*getmin(Phi[nblig-1][j+1]-Phi[nblig-1][j],0.0)+
                        getmax(Phi[nblig-1][j]-Phi[nblig-2][j],0.0)*getmax(Phi[nblig-1][j]-Phi[nblig-2][j],0.0)+
                        getmin(Phi[nblig-2][j]-Phi[nblig-1][j],0.0)*getmin(Phi[nblig-2][j]-Phi[nblig-1][j],0.0));

   }
   }
   
   else{
   for(j=1;j<nbcol-1;j++){

   norme_Phi[nblig-1][j]=sqrt(getmin(Phi[nblig-1][j]-Phi[nblig-1][j-1],0.0)*getmin(Phi[nblig-1][j]-Phi[nblig-1][j-1],0.0)+
                        getmax(Phi[nblig-1][j+1]-Phi[nblig-1][j],0.0)*getmax(Phi[nblig-1][j+1]-Phi[nblig-1][j],0.0)+
                        getmin(Phi[nblig-1][j]-Phi[nblig-2][j],0.0)*getmin(Phi[nblig-1][j]-Phi[nblig-2][j],0.0)+
                        getmax(Phi[nblig-2][j]-Phi[nblig-1][j],0.0)*getmax(Phi[nblig-2][j]-Phi[nblig-1][j],0.0));

   }

   }

   /****************************************************************/
   /*******************general case*********************************/
   /****************************************************************/
    if (regul <=0.0){
    for(i=1;i<nblig-1;i++){
      for(j=1;j<nbcol-1;j++){
        
        norme_Phi[i][j]=sqrt(getmax(Phi[i][j]-Phi[i][j-1],0.0)*getmax(Phi[i][j]-Phi[i][j-1],0.0)+
                        getmin(Phi[i][j+1]-Phi[i][j],0.0)*getmin(Phi[i][j+1]-Phi[i][j],0.0)+
                        getmax(Phi[i][j]-Phi[i-1][j],0.0)*getmax(Phi[i][j]-Phi[i-1][j],0.0)+
                        getmin(Phi[i+1][j]-Phi[i][j],0.0)*getmin(Phi[i+1][j]-Phi[i][j],0.0));


      }
    }
    }
    
    else{
    for(i=1;i<nblig-1;i++){
      for(j=1;j<nbcol-1;j++){
        
        norme_Phi[i][j]=sqrt(getmin(Phi[i][j]-Phi[i][j-1],0.0)*getmin(Phi[i][j]-Phi[i][j-1],0.0)+
                        getmax(Phi[i][j+1]-Phi[i][j],0.0)*getmax(Phi[i][j+1]-Phi[i][j],0.0)+
                        getmin(Phi[i][j]-Phi[i-1][j],0.0)*getmin(Phi[i][j]-Phi[i-1][j],0.0)+
                        getmax(Phi[i+1][j]-Phi[i][j],0.0)*getmax(Phi[i+1][j]-Phi[i][j],0.0));


      }
    }
      
    }

    /**************************************************************/
    /**************construction of the ballon force f**************/
    /**************************************************************/
    
    for(i=0;i<nblig;i++){
      for(j=0;j<nbcol;j++){
        f[i][j]=regul*dt*norme_Phi[i][j]*g[i][j];

      }
    }
    for(i=0;i<nblig;i++){
    free(norme_Phi[i]);
    }
    free(norme_Phi);
}



    /***************************************************************/
    /*************function gradient_Phi_x2**************************/
    /*****discretization of the partial derivative with respect to */
    /*****the x variable (columns).*********************************/
    /***************************************************************/

    void gradient_Phi_x2(int nblig,int nbcol,double ** Phi,double ** gradPhi_x){

    int i,j;

    /**********************right edge********************/
    for(i=0;i<nblig;i++){
    gradPhi_x[i][nbcol-1]=Phi[i][nbcol-1]-Phi[i][nbcol-2];
    }
    /**********************left edge*********************/
    for(i=0;i<nblig;i++){
    gradPhi_x[i][0]=Phi[i][1]-Phi[i][0];
    }


    /**********************general case*******************/
    for(i=0;i<nblig;i++){
      for(j=1;j<nbcol-1;j++){
      gradPhi_x[i][j]=(Phi[i][j+1]-Phi[i][j-1])/(2*h);
      }
    }



    }


    /***************************************************************/
    /*************function gradient_Phi_y2**************************/
    /*****discretization of the partial derivative with respect to*/
    /*****the y variable (rows).************************************/
    /***************************************************************/

    void gradient_Phi_y2(int nblig,int nbcol,double ** Phi,double ** gradPhi_y){

    int i,j;

    /**********************upper edge********************/
    for(j=0;j<nbcol;j++){
    gradPhi_y[0][j]=Phi[1][j]-Phi[0][j];
    }
    /**********************lower edge*********************/
    for(j=0;j<nbcol;j++){
    gradPhi_y[nblig-1][j]=Phi[nblig-1][j]-Phi[nblig-2][j];
    }


    /**********************general case*******************/
    for(i=1;i<nblig-1;i++){
      for(j=0;j<nbcol;j++){
      gradPhi_y[i][j]=(Phi[i+1][j]-Phi[i-1][j])/(2*h);
      }
    }

    }

    /*****************************************************************/
    /***************function norme_gradPhi****************************/
    /*********computes the norm of the gradient of Phi****************/
    /*****************************************************************/

    void norme_gradPhi(int nblig,int nbcol,double ** Phi,double ** norme_Phi){

    int i,j;
    double ** gradPhi_x;
    double ** gradPhi_y;
    gradPhi_x=(double**)malloc(nblig*sizeof(double*));
    gradPhi_y=(double**)malloc(nblig*sizeof(double*));
    for(i=0;i<nblig;i++){
    gradPhi_x[i]=(double*)calloc(nbcol,sizeof(double));
    gradPhi_y[i]=(double*)calloc(nbcol,sizeof(double));
    }

    gradient_Phi_x2(nblig,nbcol,Phi,gradPhi_x);
    gradient_Phi_y2(nblig,nbcol,Phi,gradPhi_y);

    for(i=0;i<nblig;i++){
      for(j=0;j<nbcol;j++){
      norme_Phi[i][j]=sqrt(gradPhi_x[i][j]*gradPhi_x[i][j]+gradPhi_y[i][j]*gradPhi_y[i][j]);
      }
    }
    
    for(i=0;i<nblig;i++){
    free(gradPhi_x[i]);
    free(gradPhi_y[i]);
    }
    free(gradPhi_x);
    free(gradPhi_y);

    }
    
    /******************************************************************/
    /***************function vectorisation_colonne*********************/
    /******concatenates the columns of the matrix given as argument****/
    /******************************************************************/

    void vectorisation_colonne(int nblig, int nbcol,double** Phi,double * v_Phi_col){
    int i,j;
    for(i=0;i<nblig;i++){
      for(j=0;j<nbcol;j++){
      v_Phi_col[j*nblig+i]=Phi[i][j];
      }
    }
    }

    /******************************************************************/
    /****************function vectorisation_ligne**********************/
    /******concatenates the rows    of the matrix given as argument****/
    /******************************************************************/

    void vectorisation_ligne(int nblig, int nbcol,double** Phi,double * v_Phi_lig){
    int i,j;
    for(i=0;i<nblig;i++){
      for(j=0;j<nbcol;j++){
      v_Phi_lig[i*nbcol+j]=Phi[i][j];
      }
    }
    }

    /********************************************************************/
    /******************function matrix_col*******************************/
    /*****creation of a matrix by filling it columns by columns**********/
    /********************************************************************/
    
    
    void matrix_col(int nblig,int nbcol,double * v_Phi,double** mat_col){
    
    int i,j;
    for(i=0;i<nblig;i++){
      for(j=0;j<nbcol;j++){
      mat_col[i][j]=v_Phi[j*nblig+i];
      }
    }

    }
    
    
    /********************************************************************/
    /******************function matrix_lig*******************************/
    /*****creation of a matrix by filling it columns by columns**********/
    /********************************************************************/
    
    
    void matrix_lig(int nblig,int nbcol,double * v_Phi,double** mat_lig){
    
    int i,j;
    for(i=0;i<nblig;i++){
      for(j=0;j<nbcol;j++){
      mat_lig[i][j]=v_Phi[i*nbcol+j];
      }
    }

    }

    /**********************************************************************/
    /*******creation of the sparse matrix of the first linear subsystem****/
    /****there is at most 3 non-zero elements on one row of the matrix*****/
    /**********************************************************************/

    void matrix_A_x2(int nblig, int nbcol,double ** g,double ** norme_Phi,double ** A_x){

    int i,j,nlig;

    /********************left edge****************/

    for(i=0;i<nblig;i++){


    nlig=i*nbcol;

    if(norme_Phi[i][0]==0||g[i][0]==0 ||g[i][1]==0){
    A_x[nlig][0]=0;
    A_x[nlig][1]=1;
    A_x[nlig][2]=0;
    }
    else{
    A_x[nlig][0]=0;
    A_x[nlig][1]=1+8*dt*norme_Phi[i][0]/((norme_Phi[i][0]/g[i][0])+(norme_Phi[i][1]/g[i][1]));
    A_x[nlig][2]=-8*dt*norme_Phi[i][0]/((norme_Phi[i][0]/g[i][0])+(norme_Phi[i][1]/g[i][1]));
    }
    }

    /******************right edge*****************/

    for(i=0;i<nblig;i++){

    nlig=i*nbcol+nbcol-1;

    if(norme_Phi[i][nbcol-1]==0||g[i][nbcol-1]==0 ||g[i][nbcol-2]==0){

    A_x[nlig][2]=0;
    A_x[nlig][1]=1;
    A_x[nlig][0]=0;
    }
    else{
    A_x[nlig][2]=0;
    A_x[nlig][1]=1+8*dt*norme_Phi[i][nbcol-1]/((norme_Phi[i][nbcol-1]/g[i][nbcol-1])+(norme_Phi[i][nbcol-2]/g[i][nbcol-2]));
    A_x[nlig][0]=-8*dt*norme_Phi[i][nbcol-1]/((norme_Phi[i][nbcol-1]/g[i][nbcol-1])+(norme_Phi[i][nbcol-2]/g[i][nbcol-2]));
    }
    }

    /***************general case*******************/
    
    for(i=0;i<nblig;i++){
      for(j=1;j<nbcol-1;j++){
      nlig=i*nbcol+j;

      if(norme_Phi[i][j]==0||g[i][j]==0||(norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i][j-1]==0 && g[i][j+1]==0)){
      A_x[nlig][0]=0;
      A_x[nlig][1]=1;
      A_x[nlig][2]=0;
      }
      else if(norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i][j-1]==0 && g[i][j+1]!=0){
      A_x[nlig][0]=0;
      A_x[nlig][1]=1+4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j+1]/g[i][j+1]));
      A_x[nlig][2]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j+1]/g[i][j+1]));
      }
      else if(norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i][j-1]!=0 && g[i][j+1]==0){
      A_x[nlig][0]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j-1]/g[i][j-1]));
      A_x[nlig][1]=1+4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j-1]/g[i][j-1]));
      A_x[nlig][2]=0;
      }
      else if(norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i][j-1]!=0 && g[i][j+1]!=0){
      A_x[nlig][0]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j-1]/g[i][j-1]));
      A_x[nlig][1]=1+4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j-1]/g[i][j-1]))+
                      4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j+1]/g[i][j+1]));
      A_x[nlig][2]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i][j+1]/g[i][j+1]));

      }
      else printf("there is an error");
      }
    }


    }


    /**********************************************************************/
    /*******creation of the sparse matrix of the second linear subsystem***/
    /****there is at most 3 non-zero elements on one row of the matrix*****/
    /**********************************************************************/

    void matrix_A_y2(int nblig, int nbcol,double ** g,double ** norme_Phi,double ** A_y){

    int i,j,nlig;

    /*******************upper edge***********************/
    for(j=0;j<nbcol;j++){
    nlig=j*nblig;

    if(norme_Phi[0][j]==0 || g[0][j]==0 || g[1][j]==0){
    A_y[nlig][0]=0;
    A_y[nlig][1]=1;
    A_y[nlig][2]=0;
    }
    else{
    A_y[nlig][0]=0;
    A_y[nlig][1]=1+8*dt*norme_Phi[0][j]/((norme_Phi[0][j]/g[0][j])+(norme_Phi[1][j]/g[1][j]));
    A_y[nlig][2]=-8*dt*norme_Phi[0][j]/((norme_Phi[0][j]/g[0][j])+(norme_Phi[1][j]/g[1][j]));
    }
    }
    
    /******************lower edge************************/
    for(j=0;j<nbcol;j++){
    nlig=j*nblig+nblig-1;

    if(norme_Phi[nblig-1][j]==0 || g[nblig-1][j]==0 || g[nblig-2][j]==0){
    A_y[nlig][0]=0;
    A_y[nlig][1]=1;
    A_y[nlig][2]=0;
    }
    else{
    A_y[nlig][2]=0;
    A_y[nlig][1]=1+8*dt*norme_Phi[nblig-1][j]/((norme_Phi[nblig-1][j]/g[nblig-1][j])+(norme_Phi[nblig-2][j]/g[nblig-2][j]));
    A_y[nlig][0]=-8*dt*norme_Phi[nblig-1][j]/((norme_Phi[nblig-1][j]/g[nblig-1][j])+(norme_Phi[nblig-2][j]/g[nblig-2][j]));
    }
    }
    /******************general case***********************/
    
    for(i=1;i<nblig-1;i++){
      for(j=0;j<nbcol;j++){

        nlig=j*nblig+i;

        if(norme_Phi[i][j]==0 ||g[i][j]==0 || (norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i+1][j]==0 && g[i-1][j]==0)){
        A_y[nlig][0]=0;
        A_y[nlig][1]=1;
        A_y[nlig][2]=0;
        }
        else if(norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i-1][j]==0 && g[i+1][j]!=0){
        A_y[nlig][0]=0;
        A_y[nlig][1]=1+4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i+1][j]/g[i+1][j]));
        A_y[nlig][2]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i+1][j]/g[i+1][j]));
        }
        else if(norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i-1][j]!=0 && g[i+1][j]==0){
        A_y[nlig][0]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i-1][j]/g[i-1][j]));
        A_y[nlig][1]=1+4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i-1][j]/g[i-1][j]));
        A_y[nlig][2]=0;
        }
        else if(norme_Phi[i][j]!=0 && g[i][j]!=0 && g[i-1][j]!=0 && g[i+1][j]!=0){
        A_y[nlig][0]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i-1][j]/g[i-1][j]));
        A_y[nlig][1]=1+4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i-1][j]/g[i-1][j]))+
                4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i+1][j]/g[i+1][j]));
        A_y[nlig][2]=-4*dt*norme_Phi[i][j]/((norme_Phi[i][j]/g[i][j])+(norme_Phi[i+1][j]/g[i+1][j]));
        }
        else printf("there is an error");
      }
    }


    }
    
    /************************Thomas Algorithm************************/
    /*The sparse matrix of the system is decomposed into the product*/
    /*LR with L a lower bidiagonal matrix with ones on the diagonal**/
    /*and R an upper bidiagonal matrix*******************************/
    /****************************************************************/

    void thomas(int nblig,double ** diag,double * second_member,double *solution){
    int k;
    double *m;
    double *l;
    double *y;

    m=(double*)malloc(nblig*sizeof(double));
    l=(double*)malloc((nblig-1)*sizeof(double));
    y=(double*)malloc(nblig*sizeof(double));


    
    m[0]=diag[0][1];

    for(k=0;k<nblig-1;k++){
    l[k]=diag[k+1][0]/m[k];
    m[k+1]=diag[k+1][1]-diag[k][2]*l[k];
    }


    /****************solution of the first linear subsystem**********/
    y[0]=second_member[0];

    for(k=1;k<nblig;k++){
    y[k]=second_member[k]-y[k-1]*l[k-1];
    }

    /****************solution of the second linear subsystem**********/

    solution[nblig-1]=y[nblig-1]/m[nblig-1];

    for(k=1;k<nblig;k++){
    solution[nblig-1-k]=(y[nblig-1-k]-diag[nblig-1-k][2]*solution[nblig-k])/m[nblig-1-k];
    }

    free(m);free(l);free(y);
    }
    
    
    
    /******************************reinitialization****************************/
    /**************************************************************************/
    /***********************by Russo and Smereka*******************************/
    /**************************************************************************/
    
    void reinitialization(int nblig, int nbcol,double **Phi0,double **Phik,double ** sgnPhi0){

    int i,j,k;
    double epsilon=0.00001;
    double gij,Dmx,Dpx,Dmy,Dpy,Dmxmax,Dmxmin,Dpxmax,Dpxmin,Dmymax,Dmymin,Dpymax,Dpymin,dij;
    double pas_temp=0.05;


    for(i=0;i<nblig;i++){
      for(j=0;j<nbcol;j++){
      sgnPhi0[i][j]=Phi0[i][j]/sqrt(Phi0[i][j]*Phi0[i][j]+epsilon);
      Phik[i][j]=Phi0[i][j];
      }
    }

    for(k=0;k<compteur;k++){


      for(i=1;i<nblig-1;i++){
        for(j=1;j<nbcol-1;j++){

             Dmx=(Phik[i][j]-Phik[i-1][j])/h;
             Dpx=(Phik[i+1][j]-Phik[i][j])/h;
             Dmy=(Phik[i][j]-Phik[i][j-1])/h;
             Dpy=(Phik[i][j+1]-Phik[i][j])/h;
             Dmxmax=getmax(Dmx,0.0);
             Dmxmin=getmin(Dmx,0.0);
             Dpxmax=getmax(Dpx,0.0);
             Dpxmin=getmin(Dpx,0.0);
             Dmymax=getmax(Dmy,0.0);
             Dmymin=getmin(Dmy,0.0);
             Dpymax=getmax(Dpy,0.0);
             Dpymin=getmin(Dpy,0.0);


             if(Phi0[i][j]>=0.0){
              gij=sqrt(getmax(Dmxmax*Dmxmax,Dpxmin*Dpxmin)+getmax(Dmymax*Dmymax,Dpymin*Dpymin))-1.0;
             }
             else if(Phi0[i][j]<0){

             gij=sqrt(getmax(Dmxmin*Dmxmin,Dpxmax*Dpxmax)+getmax(Dmymin*Dmymin,Dpymax*Dpymax))-1.0;
             }



             dij=(2*h*Phi0[i][j])/sqrt((Phi0[i+1][j]-Phi0[i-1][j])*(Phi0[i+1][j]-Phi0[i-1][j])+
                 (Phi0[i][j+1]-Phi0[i][j-1])*(Phi0[i][j+1]-Phi0[i][j-1])+epsilon);


             if((Phi0[i][j]*Phi0[i-1][j]<0.0)||(Phi0[i][j]*Phi0[i+1][j]<0.0)||(Phi0[i][j]*Phi0[i][j-1]<0.0)||(Phi0[i][j]*Phi0[i][j+1]<0.0)){

             Phik[i][j]=Phik[i][j]-(pas_temp/h)*(sgnPhi0[i][j]*getabs(Phik[i][j])-dij);

             }
             else{
             Phik[i][j]=Phik[i][j]-pas_temp*sgnPhi0[i][j]*gij;

             }


        }
      }

      /**************************boundary conditions********************/
            for(i=1;i<nblig-1;i++)
             {
             Phik[i][0]=Phik[i][1];
             Phik[i][nbcol-1]=Phik[i][nbcol-2];
             }

            for(j=1;j<nbcol-1;j++)
             {
             Phik[0][j]=Phik[1][j];
             Phik[nblig-1][j]=Phik[nblig-2][j];
             }


            Phik[0][0]=Phik[1][1];
            Phik[0][nbcol-1]=Phik[1][nbcol-2];
            Phik[nblig-1][0]=Phik[nblig-2][1];
            Phik[nblig-1][nbcol-1]=Phik[nblig-2][nbcol-2];

    }


    }
    
    
    /**********************************************************************/
    /**********C infinity regularization of the Heaviside function*********/
    /**********************************************************************/
    
    double heaviside(double eps,double phi){
    
    #define PI 3.1415926
    return(0.5*(1.0+(2/PI)*atan(phi/eps)));
    }

    /**********************************************************************/
    /**********C infinity regularization of the Dirac function*************/
    /**********************************************************************/

    double dirac(double eps,double phi){
    #define PI 3.1415926
    return((eps/PI)/(eps*eps+phi*phi));
    }
    


  /************************************************************************/
  /*********************function conv2*************************************/
  /************************************************************************/
  
  
     double conv2(int x_lig,int y_col,int level, double ** Psi,double ** contrib,double ** Phix,double ** Phiy,int window){

       int k,l;
       double prodconv1=0.0;
  
       double res;
  
       for(k=x_lig;k<x_lig+2*window+1;k++){
  
           for(l=y_col;l<y_col+2*window+1;l++){
  
           if(contrib[k][l]!=0){
  
           prodconv1=prodconv1+heaviside(epsilon,level-Psi[k][l])*heaviside(epsilon,level+Psi[k][l])*Phix[k][l]*
                     (-2.0/(1.0*(level+3)*(level+3)))*(x_lig+window-k)*exp(1.0*(-(x_lig+window-k)*(x_lig+window-k)-(y_col+window-l)*(y_col+window-l))/
                     (1.0*(level+3)*(level+3)));
  
  
           }
           }
  
       }

      res=prodconv1;
  
  
  
  
       return(res);
  
       }

  /************************************************************************/
  /*********************function conv3*************************************/
  /************************************************************************/
  
  
     double conv3(int x_lig,int y_col,int level, double ** Psi,double ** contrib,double ** Phix,double ** Phiy,int window){

     int k,l;
     double prodconv1=0.0;
     double res;


     for(k=x_lig;k<x_lig+2*window+1;k++){

         for(l=y_col;l<y_col+2*window+1;l++){

         if(contrib[k][l]!=0){

         prodconv1=prodconv1+heaviside(epsilon,level-Psi[k][l])*heaviside(epsilon,level+Psi[k][l])*Phiy[k][l]*
                   (-2.0/((level+3)*(level+3)))*(y_col+window-l)*exp(1.0*(-(x_lig+window-k)*(x_lig+window-k)-(y_col+window-l)*(y_col+window-l))/
                   (1.0*(level+3)*(level+3)));



         }
         }

     }
    res=prodconv1;




     return(res);

     }



     /*********************************************************************/
     /*********************topological force*******************************/
     /*********************************************************************/

      void topological_force(int nblig,int nbcol,double ** Phi,double ** topo,int level,int window,double ** Phix,
                             double ** Phiy,double ** Psi, double ** contrib) {

      int i,j;



      for(i=window;i<nblig+window;i++){
        for(j=window;j<nbcol+window;j++){
        Psi[i][j]=Phi[i-window][j-window];
        contrib[i][j]=1.0;
        }
      }


      /***************************************************************************/
      /********************Creation of masks containing the partial***************/
      /*************** derivative values on each node of the mesh*****************/
      /***************************************************************************/

      /********************mask Phix***************/
      /*************upper edge*******************/
      for(j=window;j<nbcol+window;j++){
      Phix[window][j]=Psi[window+1][j]-Psi[window][j];
      }
      /************lower edge********************/
      for(j=window;j<nbcol+window;j++){
      Phix[nblig+window-1][j]=Psi[nblig+window-1][j]-Psi[nblig+window-2][j];
      }
      /*************general case*****************/
      for(i=window+1;i<nblig+window-1;i++){
        for(j=window;j<nbcol+window;j++){
        Phix[i][j]=0.5*(Psi[i+1][j]-Psi[i-1][j]);
        }
      }



      /******************mask Phiy***************/
      /****************left edge*****************/
      for(i=window;i<nblig+window;i++ ){
        Phiy[i][window]=Psi[i][window+1]-Psi[i][window];
      }
      /****************right edge****************/
      for(i=window;i<nblig+window;i++ ){
        Phiy[i][window+nbcol-1]=Psi[i][window+nbcol-1]-Psi[i][window+nbcol-2];
      }
      /*****************general case*************/
      for(i=window;i<nblig+window;i++){
        for(j=window+1;j<window+nbcol-1;j++){
        Phiy[i][j]=0.5*(Psi[i][j+1]-Psi[i][j-1]);
        }
      }



      for(i=0;i<nblig;i++){
        for(j=0;j<nbcol;j++){
          if(heaviside(epsilon,level-Phi[i][j])*heaviside(epsilon,level+Phi[i][j])>0) {

        topo[i][j]=2*mu*dt*heaviside(epsilon,level-Phi[i][j])*heaviside(epsilon,level+Phi[i][j])*(conv2(i,j,level,Psi,contrib,Phix,Phiy,window)+
                  conv3(i,j,level,Psi,contrib,Phix,Phiy,window)) ;
          }

        }
      }





  }