/* Porcédure pour calculer les solutions de l'équation de Laplace 
Timothée Schmoderer
INSA Rouen Normandie
Cc 2017
*/

int Nx,Qt;

ifstream param("files/parameters");
param >> Nx >> Qt;

real[int] fdown(Nx+1);
real[int] fup(Nx+1);
real[int] mleft(Qt+1);
real[int] mright(Qt+1);
real[int,int] secondMembre(Qt+1,Nx+1); // N+1 colonnes Q+1 lignes

ifstream file0("files/fD");
ifstream file1("files/fU");
ifstream file2("files/mL");
ifstream file3("files/mR");
ifstream file4("files/d");

for (int j=0;j<Nx+1;j++) {
	file0 >> fdown(j);
	file1 >> fup(j);
}
for (int i=Qt;i>=0;i--) {
    file2 >> mleft(i);
    file3 >> mright(i);
}

for (int i=0;i<Qt+1;i++) {
    for (int j=0;j<Nx+1;j++) {
        file4 >> secondMembre(i,j);
    }
}

func fD=fdown(floor(Nx*x));
func fU=fup(floor(Nx*x));
func mL=mleft(floor(Qt*y));
func mR=mright(floor(Qt*y));
func f=secondMembre(floor(Qt*y),floor(Nx*x));

mesh Th=readmesh("files/maillage.msh"); //plot(Th,wait=1);

fespace Vh(Th,P1);
Vh uh, vh;

//		  3
//     4|	|2
//		  1

problem Poisson(uh,vh) = int2d(Th)(dx(uh)*dx(vh)+dy(uh)*dy(vh)) - int2d(Th)(f*vh) + on(4,uh=mL) + on(2,uh=mR) + on(1,uh=fD) + on(3,uh=fU);
Poisson;

plot(uh,fill=true,nbiso=50,wait=false,value=true);

ofstream output("files/solution");

for (int j=0;j<uh[].n;j++) {
	output << uh[][j] << endl;
}


