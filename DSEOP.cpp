#include <iostream>
#include <cmath>
#include "Random64.h" //Clase Random
#include <stdio.h>  //archivos
#include <fstream>  //mostrar en archivos
#include "Vector.h" //Clase Vectores
#include <omp.h>

using namespace std;
const double METRO=6100;
const int NODES=20;
const char* NAME="7.dat";
const double Tmin = 1.0, Tmax = 500.0, Tstep = 0.5;
const double Bmax = 0.3, Bstep = 0.005;
const double BBini=0.05,Tini=350.0;
const double pi = 3.14159265358979323846;
const int L=11;
const int Lz=11;
const int L2=L*L;
const int L3=L*L*Lz;
const double KU=4.1E5;
const double ratioCo=1.25E-10;
const double D0=3.544E-10+ratioCo;
const double VOL=(4.0/3.0)*pi*pow(ratioCo,3);
const double MuBohr=9.274009994E-24;
const double MuCo=3.87*MuBohr;
const double LAMBDA2=1E-8;
const double Jpd=1.393E-17;
const double wpd=0*1.7E-20;
const double Mu0=4*pi*1E-7;
const double Kb=1.3806488E-23;
const double PerCo=0.12;
const double PerO=0.54+PerCo;

vector3D DIPOLARF(int i, int j, int k, vector3D MU[L][L][Lz], int TIPO[L][L][Lz]) {
	vector3D DIPOLAR;
	DIPOLAR.cargue(0, 0, 0);
	vector3D R12;
	for (int ii=0; ii<11; ii++) {
		for (int jj=0; jj<11; jj++) {
			for (int kk=0; kk<11; kk++) {
				if (ii!=5 && jj!=5 && kk!=5&& TIPO[(i-5+ii+L)%L][(j-5+jj+L)%L][(k-5+kk+Lz)%Lz]==1) {
					R12.cargue((i-5+ii+L)%L-i, (j-5+jj+L)%L-j, (k-5+kk+Lz)%Lz-k);
					DIPOLAR+=(MU[(i-5+ii+L)%L][(j-5+jj+L)%L][(k-5+kk+Lz)%Lz]/pow(norma(R12), 3))-(3*(MU[(i-5+ii+L)%L][(j-5+jj+L)%L][(k-5+kk+Lz)%Lz]*R12)/pow(norma(R12), 5))*R12;
				}
			}
		}
	}
	//DIPOLAR.show();
	return DIPOLAR;
}

vector3D EXCHANGE(int i, int j, int k, vector3D MU[L][L][Lz], int TIPO[L][L][Lz]) {
	vector3D MUMU;
	MUMU.cargue(0, 0, 0);
	for (int ii=0; ii<3; ii++) {
		for (int jj=0; jj<3; jj++) {
			for (int kk=0; kk<3; kk++) {
				if (ii!=1 && jj!=1 && kk!=1 && TIPO[(i-1+ii+L)%L][(j-1+jj+L)%L][(k-1+kk+Lz)%Lz]==1) {
					MUMU+=MU[(i-1+ii+L)%L][(j-1+jj+L)%L][(k-1+kk+Lz)%Lz];
				}
			}
		}
	}
	//MUMU.show();
	return MUMU;
}

vector3D SUPER(int i, int j, int k, int TIPO[L][L][Lz], vector3D MU[L][L][Lz]) {
	vector3D SS;
	SS.cargue(0,0,0);
	for (int ii=0; ii<3; ii++) {
		for (int jj=0; jj<3; jj++) {
			for (int kk=0; kk<3; kk++) {
				if (ii!=1 && jj!=1 && kk!=1 && TIPO[(i-1+ii+L)%L][(j-1+jj+L)%L][(k-1+kk+Lz)%Lz]==2) {
					SS += MU[(i-1+ii+L)%L][(j-1+jj+L)%L][(k-1+kk+Lz)%Lz];
				}
			}
		}
	}
	//SS.show();
	return SS;
}

double CalculeE(double Bex, vector3D MU[L][L][Lz], vector3D e[L][L][Lz], int TIPO[L][L][Lz]) {
	double EC = 0;
	int ijkk=0;
	vector3D pos[L3];
	vector3D eijk[L3];
	vector3D MUijk[L3];
	ijkk=0;
	for (int i=0; i<L; i++) {
		for (int j=0; j<L; j++) {
			for (int k=0; k<Lz; k++) {
				if(TIPO[i][j][k]==1){
					//cout<<"cumple"<<endl;
					pos[ijkk].cargue(i,j,k);
					eijk[ijkk]=e[i][j][k];
					MUijk[ijkk]=MU[i][j][k];
					ijkk+=1;
					//cout<<ijkk<<endl;
				}
			}
		}
	}
	//cout<<ijkk<<endl;
	#pragma omp parallel for reduction(+:EC)
	for (int total=0; total<ijkk; total++) {
		EC+=-KU*VOL*(eijk[total]*MUijk[total])*(eijk[total]*MUijk[total])-MuCo*Bex*MUijk[total].z()+2*LAMBDA2*((-Jpd*MuCo*MuCo/(MuBohr*MuBohr))*(EXCHANGE(pos[total].x(), pos[total].y(), pos[total].z(), MU,TIPO)*MUijk[total])+wpd*(SUPER(pos[total].x(), pos[total].y(), pos[total].z(), TIPO, MU)*MUijk[total]))+(Mu0/(4*pi*D0*D0*D0))*MuCo*MuCo*MUijk[total]*DIPOLARF(pos[total].x(), pos[total].y(), pos[total].z(), MU,TIPO);
	}
	//cout<<EC<<endl;
	return EC;
}

double CalculeM(double Bex, vector3D MU[L][L][Lz], int TIPO[L][L][Lz]) {
	double MC = 0;
	double MUijkz[L3];
	int ijkk=0;
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			for (int k = 0; k < Lz; k++) {
				if(TIPO[i][j][k]==1){
					MUijkz[ijkk]=MU[i][j][k].z();
					ijkk++;
				}
			}
		}
	}
	#pragma omp parallel for reduction(+:MC)
	for (int total = 0; total < ijkk; total++) {
		MC += MUijkz[total];
	}
	return MuCo*MC/(L3*PerCo*VOL);
}

vector3D ROTAR(vector3D MU, Crandom & ran64) {
	vector3D MUR;
	//eje aleatorio
	int eje = (int)(3 * ran64.r());
	double alpha = 2 * pi*ran64.r();
	//roto en un eje aleatorio un angulo aleatorio
	if (eje == 0) {
		MUR.cargue(MU.x(), MU.y()*cos(alpha) - MU.z()*sin(alpha), MU.y()*sin(alpha) + MU.z()*cos(alpha));
	}
	else if (eje == 1) {
		MUR.cargue(MU.x()*cos(alpha) - MU.z()*sin(alpha), MU.y(), MU.x()*sin(alpha) + MU.z()*cos(alpha));
	}
	else {
		MUR.cargue(MU.x()*cos(alpha) - MU.y()*sin(alpha), MU.x()*sin(alpha) + MU.y()*cos(alpha), MU.z());
	}
	return MUR;
}
class SpinSystem {
private:
	vector3D MU[L][L][Lz], e[L][L][Lz];
	int TIPO[L][L][Lz];
	double E, M, Bex;
public:
	void Inicie(double Bex0, Crandom & ran64);
	double GetE(void) { return E; };
	double GetM(void) { return M; };
	void Metropolis(double T,double Bex, Crandom & ran64);
	void Muestre(void);
	void NuevoV(double Bn);
	friend double CalculeE(double Bex, vector3D MU[L][L][Lz], vector3D e[L][L][Lz], int TIPO[L][L][Lz]);
	friend double CalculeM(double Bex, vector3D MU[L][L][Lz], int TIPO[L][L][Lz]);
	friend vector3D ROTAR(vector3D MU, Crandom & ran64);
	friend vector3D EXCHANGE(int i, int j, int k, vector3D MU[L][L][Lz], int TIPO[L][L][Lz]);
	friend vector3D DIPOLARF(int i, int j, int k, vector3D MU[L][L][Lz], int TIPO[L][L][Lz]);
	friend vector3D SUPER(int i, int j, int k, double TIPO[L][L][Lz], vector3D MU[L][L][Lz]);
};
void SpinSystem::Inicie(double Bex0, Crandom & ran64) {
	double Tol = 0.0000001;
	double Norma,ranran;
	int co=0,o=0,ti=0;
	Bex = Bex0;
	for (int i = 0; i < L; i++) {
		for (int j = 0; j < L; j++) {
			for (int k = 0; k < Lz; k++) {
				ranran=ran64.r();
				//Porcentaje de dopaje
				if (ranran<PerCo) {
					//Co
					TIPO[i][j][k] = 1;
					//Eje magnetizacion facil aleatorio
					e[i][j][k].cargue(2 * ran64.r() - 1, 2 * ran64.r() - 1, 2 * ran64.r() - 1);
					//Valores aleatorios en cada dirección o alineados en z
					// /*
					MU[i][j][k].cargue(2 * ran64.r() - 1, 2 * ran64.r() - 1, 2 * ran64.r() - 1);
					//Normalizar
					Norma = norma(MU[i][j][k]);
					if (Norma > Tol) {
					MU[i][j][k] = MU[i][j][k] / Norma;
					}//*/
					//MU[i][j][k].cargue(0, 0, 1.0);
					co++;
				}
				else if (ranran<PerO) {
					//O2
					TIPO[i][j][k] = 2;
					MU[i][j][k].cargue(0, 0, 1);
					e[i][j][k].cargue(0, 0, 0);
					o++;
				}
				else {
					//another
					TIPO[i][j][k]=3;
					MU[i][j][k].cargue(0, 0, 0);
					e[i][j][k].cargue(0, 0, 0);
					ti++;
				}
				//MU[i][j][k].show();
				//cout<<norma(MU[i][j][k])<<endl;
				//cout<<TIPO[i][j][k]<<endl;
			}
		}
	}
		
	E = CalculeE(Bex, MU, e, TIPO);
	M = CalculeM(Bex, MU,TIPO);
	cout<<E<<" "<<M<<endl;
	cout<<co<<" "<<o<<" "<<ti<<endl;
}
//L3PASOS
void SpinSystem::Metropolis(double T, double Bex, Crandom & ran64) {
	int n, i, j, k,c1=0,c2=0,c3=0;
	double dE, Eold;
	vector3D MUN, MUO;
	for (int mcs=0; mcs<L3; mcs++) { //Un MCSS
		//  Escojo un dipolo al azar;
		n=(int)L2*ran64.r(); i=n/L; j=n%L; k=(int)Lz*ran64.r();
		//Calculo la nueva energia que se produciría si lo volteo con una matriz aleatoria;
		MUO=MU[i][j][k];
		if (TIPO[i][j][k]==1) {
			Eold = E;
			//cout<<E<<endl;
			MUN = ROTAR(MU[i][j][k], ran64);
			MU[i][j][k] = MUN;
			dE = CalculeE(Bex, MU, e, TIPO) - Eold;
			//cout<<exp(-dE/(Kb*T))<<endl;
			//cout<<dE<<" "<<Eold<<endl;
			if (dE<=0.0) {
				E += dE; M = CalculeM(Bex, MU,TIPO);
				c1++;
				//cout<<"c1 "<<i<<" "<<j<<" "<<k<<endl;
			}//lo volteo;
			else if (ran64.r()<exp(-METRO*dE/(Kb*T))) {
				E += dE; M = CalculeM(Bex, MU,TIPO);
				c2++;
				//cout<<"c2 "<<i<<" "<<j<<" "<<k<<endl;
				//cout<<exp(-dE/(Kb*T))<<" "<<dE/(Kb*T)<<endl;
				//cout<<dE<<" "<<Eold<<" "<<E<<endl;
			}//lo volteo;
			else {
				MU[i][j][k] = MUO;
				c3++;
				//cout<<"c3"<<endl; 
			}//no lo volteo
		}
	}
//	cout<<c1<<" "<<c2<<" "<<c3<<endl;
}
void SpinSystem::Muestre(void) {
	for (int i = 0; i<L; i++) {
		for (int j = 0; j < L; j++) {
			for (int k = 0; k < Lz; k++) {
				MU[i][j][k].show();
			}
		}
	}
}
void SpinSystem::NuevoV(double Bn) {
	Bex = Bn;
	E = CalculeE(Bex, MU, e, TIPO);
	M = CalculeM(Bex, MU, TIPO);
}
//equilibrium constants
const int teq = 500;// (int)(600 * pow(L3 / 8.0, 2.125));
const int tcorr = teq / 10;
const int Nmuestras = 10;
//
//MAIN PROGRAM
int main(void) {
	omp_set_num_threads(NODES);
	ofstream myfile(NAME);
	SpinSystem Ising;
	Crandom ran2(546);
	double Mprom, M, E,T,BB;
myfile<<"NODES "<<NODES<<endl;
myfile<<"First-Temperature "<<Tmin<<endl;
myfile<<"Last-Temperature "<<Tmax<<endl;
myfile<<"Temperature-Step "<<Tstep<<endl;
myfile<<"Max-External-Field "<<Bmax<<endl;
myfile<<"External-Field-Step "<<Bstep<<endl;
myfile<<"External-Field-ZFCFC "<<BBini<<endl;
myfile<<"Temperature-for-HM-curves "<<Tini<<endl;
myfile<<"L "<<L<<endl;
myfile<<"Lz "<<Lz<<endl;
myfile<<"Ku "<<KU<<endl;
myfile<<"Ratio-of-Co "<<ratioCo<<endl;
myfile<<"Interatomic-Distance "<<D0<<endl;
myfile<<"Volume-of-Co "<<VOL<<endl;
myfile<<"Bohr-Magneton "<<MuBohr<<endl;
myfile<<"Co-Magnetization "<<MuCo<<endl;
myfile<<"Lambda2 "<<LAMBDA2<<endl;
myfile<<"Jpd "<<Jpd<<endl;
myfile<<"Wpd "<<wpd<<endl;
myfile<<"Mu0 "<<Mu0<<endl;
myfile<<"Kb "<<Kb<<endl;
myfile<<"Co-Percentage "<<PerCo<<endl;
myfile<<"O-Percentage "<<PerO-PerCo<<endl;
	/*
	//PARA ESTUDIAR EL TIEMPO DE EQUILIBRIO
	myfile<<"T "<<Tini<<" BB "<<BBini<<endl;
	//Inicio el sistema
	Ising.Inicie(BBini, ran2);
	cout<<"INICIO"<<endl;
	double EE=0, EEN=0, Es2=0, Es2N=0;
	double MM=0, MMN=0, Ms2=0, Ms2N=0;
	for (int j=0; j<teq*10000; j++){
		for (int i=0; i<Nmuestras; i++) {
			Ising.Metropolis(Tini, BBini, ran2);
			EE+=Ising.GetE();
			MM+=Ising.GetM();
			Es2+=(Ising.GetE()-EEN)*(Ising.GetE()-EEN);
			Ms2+=(Ising.GetM()-MMN)*(Ising.GetM()-MMN);
		}
		EEN=EE/Nmuestras;
		EE=0;
		Es2N=Es2/Nmuestras;
		Es2=0;
		MMN=MM/Nmuestras;
		MM=0;
		Ms2N=Ms2/Nmuestras;
		Ms2=0;
		cout<<"Eprom "<<EEN<<" Edes "<<sqrt(Es2N)<<" Mprom "<<MMN<<" Mdes "<<sqrt(Ms2N)<<endl;		
		myfile<<"Eprom "<<EEN<<" Edes "<<sqrt(Es2N)<<" Mprom "<<MMN<<" Mdes "<<sqrt(Ms2N)<<endl;
	}
	 // */
	// /*
	//REALIZAR CURVAS DE HISTERESIS
	//Inicio el sistema
	Ising.Inicie(0, ran2);
	for (BB=0; BB<Bmax; BB+=Bstep) {
		Ising.NuevoV(BB);
		cout<<"BB  "<<BB<<endl;
		//Equilibro
		for (int t=0; t<teq; t++) {
			//cout<<t<<endl;
			Ising.Metropolis(Tini, BB, ran2);
		}
		cout<<"eq  "<<endl;
		//Inicio Acumuladores en cero
		Mprom=0, M=0, E=0;
		//Tomo muestras
		for (int t=0; t<Nmuestras; t++) {
			M=Ising.GetM();
			E+=Ising.GetE();
			Mprom+=M;
			for (int ccorr=0; ccorr<tcorr; ccorr++) {
				Ising.Metropolis(Tini, BB, ran2);
			}
		}
		cout<<"Finish"<<endl;
		//Normalizo los acumuladores
		Mprom/=Nmuestras;
		//Imprimo
		cout<<Mprom<<" "<<E/Nmuestras<<endl;
		myfile<<Tini<<" "<<BB<<" "<<Mprom<<endl;
	}
	for (BB=Bmax+0.5*Bstep; BB>-Bmax; BB-=Bstep) {
		Ising.NuevoV(BB);
		cout<<"BB  "<<BB<<endl;
		//Equilibro
		for (int t=0; t<teq; t++) {
			Ising.Metropolis(Tini, BB, ran2);
		}
		cout << "eq  " << endl;
		//Inicio Acumuladores en cero
		Mprom = 0, M = 0, E = 0;
		//Tomo muestras
		for (int t=0; t<Nmuestras; t++) {
			M=Ising.GetM();
			Mprom+=M;
			for (int ccorr = 0; ccorr < tcorr; ccorr++) {
				Ising.Metropolis(Tini, BB, ran2);
			}
		}
		cout << "Finish" << endl;
		//Normalizo los acumuladores
		Mprom/=Nmuestras;
		//Imprimo
		myfile<<Tini<<" "<<BB<<" "<<Mprom<<endl;
		cout<<Mprom<<endl;
	}
	for (BB=-Bmax-0.6*Bstep; BB<Bmax; BB+=Bstep) {
		Ising.NuevoV(BB);
		cout << "BB  " << BB << endl;
		//Equilibro
		for (int t = 0; t<teq; t++) {
			Ising.Metropolis(Tini, BB, ran2);
		}
		cout<<"eq  "<<endl;
		//Inicio Acumuladores en cero
		Mprom=0, M=0, E=0;
		//Tomo muestras
		for (int t=0; t<Nmuestras; t++) {
			M=Ising.GetM();
			Mprom+=M;
			for (int ccorr=0; ccorr<tcorr; ccorr++) {
				Ising.Metropolis(Tini, BB, ran2);
			}
		}
		cout<<"Finish"<<endl;
		//Normalizo los acumuladores
		Mprom/=Nmuestras;
		//Imprimo
		myfile<<Tini<<" "<<BB<<" "<<Mprom<<endl;
	}
	// */
	
	 /*
	Ising.Inicie(0, ran2);
	//CERO FIEL COOL FIEL COOL

	for (int t = 0; t<teq*10; t++) {
		Ising.Metropolis(Tmin, 0, ran2);
	}
	cout << "Feq  " << endl;
	BB = BBini;
	Ising.NuevoV(BB);
	for (int t=0; t<teq*10; t++) {
		Ising.Metropolis(Tmin, BB, ran2);
	}
	cout<<"Seq  "<<endl;
	for (T=Tmin; T<Tmax; T+=Tstep) {
		cout << "T  " << T << endl;
		//Equilibro
		for (int t=0; t<teq; t++) {
			Ising.Metropolis(T, BB, ran2);
		}
		cout << "eq  " << endl;
		//Inicio Acumuladores en cero
		Mprom=0, M=0, E=0;
		//Tomo muestras
		for (int t=0; t<Nmuestras; t++) {
			M = Ising.GetM();
			Mprom+=M;
			for (int ccorr = 0; ccorr < tcorr; ccorr++) {
				Ising.Metropolis(T, BB, ran2);
			}
		}
		cout << "Finish" << endl;
		//Normalizo los acumuladores
		Mprom/=Nmuestras;
		//Imprimo
		myfile << T << " " << BB << " " << Mprom << endl;
	}
	for (T=Tmax+0.5*Tstep; T>Tmin; T-=Tstep){
		cout << "T  " << T << endl;
		//Equilibro
		for (int t = 0; t<teq; t++){
			Ising.Metropolis(T, BB, ran2);
		}
		cout << "eq  " << endl;
		//Inicio Acumuladores en cero
		Mprom = 0, M = 0, E = 0;
		//Tomo muestras
		for (int t = 0; t<Nmuestras; t++) {
			M = Ising.GetM();
			Mprom += M;
			for (int ccorr = 0; ccorr < tcorr; ccorr++){
				Ising.Metropolis(T, BB, ran2);
			}
		}
		cout << "Finish" << endl;
		//Normalizo los acumuladores
		Mprom /= Nmuestras;
		//Imprimo
		myfile << T << " "  << BB << " " << Mprom << endl;
	}
	// */

	myfile.close();
	return 0;
}
