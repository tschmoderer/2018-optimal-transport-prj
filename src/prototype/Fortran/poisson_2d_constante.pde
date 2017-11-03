/*
Procédure de calcul de l'équation de Poisson dans le carré [0,1]x[0,1]
Timothée Schmoderer
INSA Rouen Normandie
Cc 2017
*/

int Nx,Qt; // Taille de la grille 

ifstream param("files/parameters");
param >> Nx >> Qt;

// les conditions aux fontières 

real[int] finit(Nx+1); // condition initial --> ira en bas
real[int] ffinal(Nx+1); // condition terminal --> ira en haut 
ifstream file0("files/f0");
ifstream file1("files/f1");
for (int j=0;j<Nx+1;j++) {
	file0 >> finit(j);
	file1 >> ffinal(j);
}

func f0=finit(floor(Nx*x));
func f1=ffinal(floor(Nx*x));

mesh Th=square(Nx,Qt); 
// plot(Th,wait=1);

fespace Vh(Th,P1);
Vh uh, vh;

//		  3
//     4|	|2
//		  1  


problem Poisson(uh,vh) = int2d(Th)(dx(uh)*dx(vh)+dy(uh)*dy(vh)) + on(1,uh=f0) + on(2,4,uh=0) + on(3,uh=f1);
Poisson;

// plot(uh,wait=1,nbiso=30,fill=true,value=true);

ofstream output("files/Y");

for (int j=0;j<uh[].n;j++) {
	output << uh[][j] << endl;
}

savemesh(Th,"files/maillage.msh");
