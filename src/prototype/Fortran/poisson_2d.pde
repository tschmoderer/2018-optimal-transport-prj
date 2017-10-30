// https://arxiv.org/pdf/1205.1293.pdf
// manuel freefem p244

/*

*/

int Nx,Qt;
int [int] d(2);

ifstream param("files/parameters");
param >> Nx >> Qt;

real[int] fdown(Nx+1);
real[int] fup(Nx+1);
real[int] mleft(Qt+1);
real[int] mright(Qt+1);
real[int,int] secondMembre(Nx+1,Qt+1);

ifstream file0("files/f0");
ifstream file1("files/f1");
ifstream file2("files/mL");
ifstream file3("files/mR");
ifstream file4("files/d");

for (int j=0;j<Nx+1;j++) {
	file0 >> fdown(j);
	file1 >> fup(j);
}
for (int i=0;i<Qt+1;i++) {
    file2 >> mleft(i);
    file3 >> mright(i);
}
for (int i=0;i<Nx+1;i++) {
    for (int j=0;j<Qt+1;j++) {
        file4 >> secondMembre(i,j);
    }
}

func fD=fdown(Nx*x);
func fU=fup(Nx*x);
func mL=mleft(Qt*x);
func mR=mright(Qt*x);
func f=secondMembre(Nx*x,Qt*y);

mesh Th=readmesh("files/maillage.msh");
//plot(Th,wait=1);

fespace Vh(Th,P1);

Vh uh, vh;

//		  3
//      4|	|2
//		  1
//


problem Poisson(uh,vh) = int2d(Th)(dx(uh)*dx(vh)+dy(uh)*dy(vh)) - int2d(Th)(f*vh) + on(2,uh=mL) + on(4,uh=mR) + on(1,uh=fD) + on(3,uh=fU);

Poisson;

plot(uh,wait=1,fill=true,value=true);

ofstream output("files/solution");
for(int j=0; j<uh[].n;j++) {
    output << uh[][j] << endl;
}
