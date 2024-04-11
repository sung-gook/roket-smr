#include <iostream>
#include <fstream>
#include<iomanip>
#include <cmath>

using namespace std;
double Ve=100., Me=0.1, M0=1.0, g = 9.81;

double dVrdt(double t, double Vr) {
	return Me *Ve / (M0 - Me * t) - g;
}

int main() {


int nStep, nPrecision=5, nWidth=15;
double *Vrk, *Veu, *VTh;
double time, dTime=0.02, timeEnd=1.;
double errorEU, errorRK;

ofstream fout;
fout.open("Ex2-3.out");

nStep = int (timeEnd/dTime) + 1;

Vrk = new double[nStep];
Veu = new double[nStep];
VTh = new double[nStep];

time = 0.0;
Vrk[0] = 0.;
Veu[0] = 0.;
VTh[0] = 0.;

cout.precision(nPrecision);
fout.precision(nPrecision);


cout << setw(nWidth) << "time" << setw(nWidth) << "V(Exact)\t" << setw(nWidth) << "V(Euler)\t" << setw(nWidth)  << "V(Runge-Kutta)\t" << setw(nWidth) << "error(Euler)\t" << "error(R-K)" << endl;
fout << setw(nWidth) << "time" << setw(nWidth) << "V(Exact)\t" << setw(nWidth) << "V(Euler)\t" <<  setw(nWidth) << "V(Runge-Kutta)\t" << setw(nWidth) << "error(Euler)\t" << "error(R-K)" << endl;
double k1, k2;
for (int i=0;i<nStep-1;i++) {
	//오일러법
	Veu[i + 1] = Veu[i] + dTime * dVrdt(time, Veu[i]);
	//룬지쿠타
	k1 = dTime * dVrdt(time, Vrk[i]);
	k2 = dTime * dVrdt(time + dTime, Vrk[i] + k1);
	time +=dTime;
	Vrk[i + 1] = Vrk[i] + 0.5*(k1 + k2);
	//이론해
	VTh[i + 1] = Ve * log(M0 / (M0 - Me * time)) - g * time;
	//오차
	errorEU = fabs((Veu[i + 1] - VTh[i + 1]) / VTh[i + 1])*100;
	errorRK = fabs((Vrk[i + 1] - VTh[i + 1]) / VTh[i + 1])*100;


	cout << setw(nWidth) << time << setw(nWidth) << VTh[i + 1] << setw(nWidth) << Veu[i + 1] << setw(nWidth) << Vrk[i + 1] << setw(nWidth)<< errorEU<< setw(nWidth) <<errorRK<< endl;
	fout << setw(nWidth) << time << setw(nWidth) << VTh[i + 1] << setw(nWidth) << Veu[i + 1] << setw(nWidth) << Vrk[i + 1] << setw(nWidth)<< errorEU<< setw(nWidth) <<errorRK<< endl;
	
  }

fout.close();
return 0;
}


