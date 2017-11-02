/*
Calcul de la constante dans la solution
genere le maillage
*/

int Nx,Qt; // Taille de la grille 

ifstream param("files/parameters");
param >> Nx >> Qt;

// les conditions aux fontières 

real[int] finit(Nx+1);
real[int] ffinal(Nx+1);
ifstream file0("files/f0");
ifstream file1("files/f1");
for (int j=0;j<Nx+1;j++) {
	file0 >> finit(j);
	file1 >> ffinal(j);
}

func f0=finit(Nx*x);
func f1=ffinal(Nx*x);

mesh Th=square(Nx,Qt); //plot(Th,wait=1);

fespace Vh(Th,P1);
Vh uh, vh;

//		  3
//     4|	|2
//		  1  


problem Poisson(uh,vh) = int2d(Th)(dx(uh)*dx(vh)+dy(uh)*dy(vh)) + on(2,4,uh=0) + on(1,uh=f0) + on(3,uh=f1);
Poisson;

plot(uh,wait=1,fill=true,value=true);

ofstream output("files/Y");
// Affichage à l'envers pour lecture fortran dans le bon sens
for(int j=uh[].n-1; j>=0;j--) {
    output << uh[][j] << endl; // a voir s'il faudrait pas tout ecrire sur une ligne 
}

savemesh(Th,"files/maillage.msh");