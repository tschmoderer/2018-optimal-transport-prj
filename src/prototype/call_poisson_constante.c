#include "mex.h" 
#include <stdio.h>
#include <string.h>

void mexFunction(int nlhs, mxArray *plhs[],	    /* Sorties */
                 int nrhs, const mxArray *prhs[])   /* Entrées */
{
/*
    nlhs : Nombre de sorties demandées
    plhs : Sorties	
                plhs[0] : sortie 1 
                plhs[1] : sortie 2 
                ...
                plhs[nlhs-1] : sortie nlhs 

    nrhs : Nombre d'entrées
    prhs : Entrées	
                prhs[0] : entrée 1 
                prhs[1] : entrée 2 
                ...
                prhs[nrhs-1] : entrée nrhs 
*/

    /* code de la fonction */
     mexPrintf("Hello World !\n");
   char command[50];

   strcpy(command, "FreeFem++ poisson_2d.pde");
   int i = system(command);
   
   nlhs = 1;
   plhs[0] = 0;

}