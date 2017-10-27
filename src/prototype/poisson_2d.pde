// https://arxiv.org/pdf/1205.1293.pdf
// manuel freefem p244

/*
Calcul de la constante dans la solution
genere le maillage
*/
load "lapack"

int Nx,Qt;
int [int] d(2);
real mu = 0.1;
real s = 0.005;

ifstream param("files/parameters.txt");
param >> Nx >> Qt;

//func f0=exp(-0.5*((x-0.1)/s)^2)/(s*sqrt(2*pi));
//func f1=exp(-0.5*((x-0.9)/s)^2)/(s*sqrt(2*pi));

real[int] finit(Nx+1);
real[int] ffinal(Nx+1);
ifstream file0("files/f0.txt");
ifstream file1("files/f1.txt");
for (int j=0;j<Nx+1;j++) {
	file0 >> finit(j);
	file1 >> ffinal(j);
}

func f0=finit(Nx*x);
func f1=ffinal(Nx*x);


mesh Th=square(Nx,Qt);
plot(Th,wait=1);

fespace Vh(Th,P1);

Vh uh, vh;

//		  3
//      4|	|2
//		  1
//


problem Poisson(uh,vh) = int2d(Th)(dx(uh)*dx(vh)+dy(uh)*dy(vh)) + on(2,4,uh=0) + on(1,uh=f0) + on(3,uh=f1);

Poisson;

plot(uh,wait=1,fill=true,value=true);

ofstream output("files/Y.txt");
for(int j=0; j<uh[].n;j++) {
    output << uh[][j]<<endl;
}

savemesh(Th,"files/maillage.msh");

