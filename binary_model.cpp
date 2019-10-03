#include <iostream>
#include <vector>
#include <string>
#include <deque>
#include <set>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <sstream>

using namespace std;
double const pi = 3.14159265;
double const e = 2.71828182;

double const T = 1800*1000.0;
double const h = 0.01;

int const NE = 2500;
int const NI = 500;
int const N = NE + NI;

double const cEE = 0.2;
double const cIE = 0.2;
double const cEI = 0.5;
double const cII = 0.5;

double JEE = 0.15;
double JEEh = 0.15;
//double JEI = 0.22;
//double JEIo = 0.15;
double JIE = 0.15;
double JII = 0.06;
double sigJ = 0.3;

double Jtmax = 0.25; //0.25
double Jtmin = 0.01;

double Jmax = 5.0*JEE;
double Jmin = 0.01*JEE;

double const hE = 1.0;
double hI = 1.0;

double const IEex = 2.0;
double const IIex = 0.5;
double mex = 0.3;
double sigex = 0.1;

double const tmE = 5.0;
double const tmI = 2.5;

int const SNE = (int)floor(NE*h/tmE + 0.001);
int const SNI = (int)floor(NI*h/tmI + 0.001);

//Short-Term Depression
double trec = 600.0;
//double usyn = 0.1;
double Jepsilon = 0.001;

//STDP
double const tpp = 20.0;
double const tpd = 40.0;
double const twnd = 600.0;
double const Cp = 0.10*JEE; //0.1*JEE;
double const Cd = Cp*tpp/tpd;
double const g = 1.25;

//homeostatic
//double hsig = 0.001*JEE/sqrt(10.0);
double hsig = 0.001*JEE;
int itauh = 100;

double hsd = 0.1;
double hh = 10.0;

double Ip = 1.0; //1.0;
//double a = 0.20;


double o1th = 0.01;

//initial input
double xEinit = 0.02;
double xIinit = 0.01;
double tinit = 100.0;

double tdur = 1000.0;

vector<double> dvec;
vector<int> ivec;
deque<int> ideque;

double dice(){
	return rand()/(RAND_MAX + 1.0);
}

double ngn(){
	double u = dice(); double v = dice();
	return sqrt(-2.0*log(u))*cos(2.0*pi*v);
}

vector<int> rnd_sample(int ktmp, int Ntmp){ // when ktmp << Ntmp
	vector<int> smpld; int xtmp; bool tof;
	while( smpld.size() < ktmp ){
		xtmp = (int)floor( Ntmp*dice() ); tof = true;
		for(int i = 0; i < smpld.size(); i++ ){
			if( xtmp == smpld[i] ) tof = false;
		}
		if(tof) smpld.push_back(xtmp);
	}
	return smpld;
}

double fd(double x, double alpha){
	return log(1.0 + alpha*x)/log(1.0 + alpha);
}

void calc(vector< vector<double> > Jo, double alpha, double usd, double JEI, double a, double Jp, int p, int ik){
	int ialpha = (int)floor(alpha + 0.01);
	int iusd = (int)floor(usd*100.1);
	int iJEI = (int)floor(JEI*1000.01);
	int ia = (int)floor(a*100.01);
	int iJp = (int)floor(Jp*100.01);

	//ostringstream ossr; //spiking data
	//ossr << "binary_model_sr_al" << ialpha <<"_u"<< iusd <<"_i"<< iJEI <<"_a"<< ia <<"_j"<< iJp <<"_p"<< p <<"_k"<< ik << ".txt"; 
	//string fstrr = ossr.str(); ofstream ofsr; ofsr.open( fstrr.c_str() );
	//ofsr.precision(10);
	//ostringstream ossw; //synaptic weight matrix
	//ossw << "binary_model_sw_al" << ialpha <<"_u"<< iusd <<"_i"<< iJEI <<"_a"<< ia <<"_j"<< iJp <<"_p"<< p <<"_k"<< ik << ".txt"; 
	//string fstrw = ossw.str(); ofstream ofsw; ofsw.open( fstrw.c_str() );
	ostringstream ossd; //mean synaptic weights among assemblies
	ossd << "binary_model_sd_al" << ialpha <<"_u"<< iusd <<"_i"<< iJEI <<"_a"<< ia <<"_j"<< iJp <<"_p"<< p <<"_k"<< ik << ".txt"; 
	string fstrd = ossd.str(); ofstream ofsd; ofsd.open( fstrd.c_str() );
	ostringstream osspf; //population firing rate
	osspf << "binary_model_spf_al" << ialpha <<"_u"<< iusd <<"_i"<< iJEI <<"_a"<< ia <<"_j"<< iJp <<"_p"<< p <<"_k"<< ik << ".txt"; 
	string fstrpf = osspf.str(); ofstream ofspf; ofspf.open( fstrpf.c_str() );

	double tauh = itauh*1000.0; int NEa = (int)floor(NE*a+0.01);

	vector<int> ptn_inv;
	for(int q = 0; q < p; q++){
		for(int i = q*NEa; i < (q+1)*NEa; i++) ptn_inv.push_back(q+1);
	}
	for(int i = p*NEa; i < NE; i++) ptn_inv.push_back(0);
	vector< vector<double> > wqqcnt;
	for(int i = 0; i < p+1; i++){
		wqqcnt.push_back(dvec);
		for(int i2 = 0; i2 < p+1; i2++) wqqcnt[i].push_back(0.0);
	}

	vector<double> ys;
	for(int i = 0; i < NE; i++) ys.push_back( 1.0/(1.0 + usd*0.05*trec/tmE) );

	vector< vector<int> > Jinidx;
	for(int i = 0; i < NE; i++){
		Jinidx.push_back(ivec);
		for(int j = 0; j < NE; j++){
			if( Jo[i][j] > Jepsilon ){
				Jinidx[i].push_back( j );
				wqqcnt[ ptn_inv[i] ][ ptn_inv[j] ] += 1.0;
			}
		}
	}

	vector<int> x;
	for(int i = 0; i < N; i++) x.push_back(0);
	set<int> spts;
	for(int i = 0; i < N; i++){
		if( i < NE && dice() < xEinit ){
			spts.insert(i); x[i] = 1;
		}
		if( i >= NE && dice() < xIinit ){
			spts.insert(i); x[i] = 1;
		}
	}

	vector<double> pfrs; double dpfr = 1000.0/(10.0*NEa);
	for(int q = 0; q < p+2; q++){
		pfrs.push_back(0.0);
	}

	vector< deque<int> > dspts;
	for(int i = 0; i < NE; i++) dspts.push_back( ideque );

	int tidx = -1; bool trtof = false;
	double u; int j; vector<int> smpld;
	set<int>::iterator it;
	double k1,k2,k3,k4; bool Iptof = true;
	for(double t = 0; t < T+h; t += h){
		smpld = rnd_sample(SNE,NE);
		for(int iidx = 0; iidx < smpld.size(); iidx++){
			int i = smpld[iidx];
			if( x[i] == 1 ){
				ys[i] -= usd*ys[i];
				it = spts.find( i );
				if( it != spts.end() ) spts.erase( it++ );
				x[i] = 0;
			}
			u = -hE + IEex*(mex + sigex*ngn());
			/* 
			if( i < 2*NEa ){
				for(int cidx = 0; cidx < cstim; cidx++){
					if( tcstims[cidx] <= t && t < tcstims[cidx] + 100.0 ){
						u += 0.5*Ip;					
					}
				}
			}
			*/
			it = spts.begin();
			while( it != spts.end() ){
				if( *it < NE){
					u += ys[*it]*Jo[i][ *it ];
				}else{
					u += Jo[i][ *it ];
				}
				++it;
			}
			if( u > 0 ){
				spts.insert(i); dspts[i].push_back(t); x[i] = 1; 
				//if( trtof ) ofsr << t << " " << i << endl;
				pfrs[i/NEa] += dpfr;
				//E-pre
				for(int ip = 0; ip < NE; ip++){
					if( Jo[ip][i] > Jepsilon && t > tinit ){
						for(int sidx = 0; sidx < dspts[ip].size(); sidx++){
							Jo[ip][i] += Cp*g*exp( -(t-dspts[ip][sidx])/tpp );
							Jo[ip][i] -= Cd*fd(Jo[ip][i]/JEE,alpha)*exp( -(t-dspts[ip][sidx])/tpd );
						}
						if( Jo[ip][i] < Jmin ) Jo[ip][i] = Jmin;
						if( Jo[ip][i] > Jmax ) Jo[ip][i] = Jmax;
					}
				}
				//E-post
				for(int jidx = 0; jidx < Jinidx[i].size(); jidx++){
					j = Jinidx[i][jidx];
					if( t > tinit){
						for(int sidx = 0; sidx < dspts[j].size(); sidx++){
							Jo[i][j] += g*Cp*exp( -(t-dspts[j][sidx])/tpp );
							Jo[i][j] -= Cd*fd(Jo[i][j]/JEE,alpha)*exp( -(t-dspts[j][sidx])/tpd );
						}
						if( Jo[i][j] < Jmin ) Jo[i][j] = Jmin;
						if( Jo[i][j] > Jmax ) Jo[i][j] = Jmax;
					}
				}
			}
		}
	
		smpld = rnd_sample(SNI,NI);
		for(int iidx = 0; iidx < smpld.size(); iidx++){
			int i = NE + smpld[iidx];
			if( x[i] == 1 ){
				it = spts.find( i );
				if( it != spts.end() ) spts.erase( it++ );
				x[i] = 0;
			}
			u = -hI + IIex*(mex + sigex*ngn()); 
			it = spts.begin();
			while( it != spts.end() ){
				u += Jo[i][ *it ]; ++it;
			}
			if( u > 0 ){
				spts.insert(i); x[i] = 1; 
				//if( trtof ) ofsr << t << " " << i << endl;
			}
		}
		if( ( (int)floor(t/h) )%10 == 0 ){
			//STD
			for(int i = 0; i < NE; i++){
				k1 = (1.0 - ys[i])/trec; k2 = (1.0 - (ys[i]+0.5*hsd*k1))/trec;
				k3 = (1.0 - (ys[i]+0.5*hsd*k2))/trec; k4 = (1.0 - (ys[i]+hsd*k3))/trec;
				ys[i] += hsd*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
			}
		}
		if( ( (int)floor(t/h) )%1000 == 0 ){
			//Homeostatic Depression
			for(int i = 0; i < NE; i++){
				for(int jidx = 0; jidx < Jinidx[i].size(); jidx++){
					j = Jinidx[i][jidx];
					k1 = (JEEh - Jo[i][j])/tauh; k2 = (JEEh - (Jo[i][j]+0.5*hh*k1))/tauh;
					k3 = (JEEh - (Jo[i][j] + 0.5*hh*k2))/tauh; k4 = (JEEh - (Jo[i][j] + hh*k3))/tauh;
					Jo[i][j] += hh*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0 + hsig*ngn();
					if( Jo[i][j] < Jmin ) Jo[i][j] = Jmin;
					if( Jo[i][j] > Jmax ) Jo[i][j] = Jmax;
				}
			}
			//boundary condition
			for(int i = 0; i < NE; i++){
				double Jav = 0.0;
				for(int jidx = 0; jidx < Jinidx[i].size(); jidx++) Jav += Jo[i][ Jinidx[i][jidx] ];
				Jav = Jav/( (double)Jinidx[i].size() );
				if( Jav > Jtmax ){
					for(int jidx = 0; jidx < Jinidx[i].size(); jidx++){
						j = Jinidx[i][jidx]; Jo[i][j] -= (Jav-Jtmax);
						if( Jo[i][j] < Jmin ) Jo[i][j] = Jmin;
					}
				}else if( Jav < Jtmin ){
					for(int jidx = 0; jidx < Jinidx[i].size(); jidx++){
						j = Jinidx[i][jidx]; Jo[i][j] += (Jtmin-Jav);
						if( Jo[i][j] > Jmax ) Jo[i][j] = Jmax;
					}
				}
			}
			//STDP
			for(int i = 0; i < NE; i++){
				for(int sidx = 0; sidx < dspts[i].size(); sidx++){
					if( t - dspts[i][0] > twnd ) dspts[i].pop_front();
				}
			}

			//population firing rate
			for(int q = 0; q < p+2; q++){
				ofspf << pfrs[q] << " "; pfrs[q] = 0.0;
			}
			ofspf << endl;
		}

		if( ( (int)floor(t/h) )%(1000*100) == 0 ){
			tidx++;
			//if( tidx < 50 || tidx > 1750){ //CAUTION!!! HEAVY DATA!!!
			//if( tidx < 5 || tidx > 1795 ){
			//	trtof = true;
			//}else{
			//	trtof = false;
			//}
			if( ( (int)floor(t/h) )%(1000*100) == 0 ){
				vector< vector<double> > wqq;
				for(int i = 0; i < p+1; i++){
					wqq.push_back( dvec );
					for( j = 0; j < p+1; j++) wqq[i].push_back(0.0);
				}
				for(int i = 0; i < NE; i++){
					for(j = 0; j < NE; j++){
						if( abs(Jo[i][j]) > Jepsilon ) wqq[ ptn_inv[i] ][ ptn_inv[j] ] += Jo[i][j];
					}
				}
				for(int i = 0; i < p+1; i++){
					for( j = 0; j < p+1; j++) ofsd << wqq[i][j]/wqqcnt[i][j] << " ";
				}
				ofsd << endl;
			}
			/*
			if( tidx == 0 || tidx == 900 || tidx == 1799 ){
				for(int i = 0; i < N; i++){
					for(j = 0; j < N; j++){
						if( abs(Jo[i][j]) > Jepsilon ) ofsw << Jo[i][j] << " " << j << " ";
					}
					ofsw << endl;
				}
			}
			*/
			
			if( ( (int)floor(t/h) )%(10000*100) == 0 ) cout << t/1000.0 << endl;
			int s = 0; it = spts.begin();
			while( it != spts.end() ){
				++s; ++it;
			}
			//cout << s << endl;
			if( s == 0 || (s > 1.0*NE && t > 200.0) ) break;
		}	
	}

}

//Patterned initialization of the weight matrix J
vector< vector<double> > calc_J(double Jp, double Jb, double JEI, double a, int p){
	int NEa = (int)floor(NE*a+0.01);

	vector< vector<double> > J;	
	int mcount = 0; 
	for(int i = 0; i < NE; i++){
		J.push_back(dvec);
		for(int j = 0; j < NE; j++){
			J[i].push_back(0.0);
			if( i != j && dice() < cEE ){
				J[i][j] += Jb*(1.0 + sigJ*ngn());
				if( J[i][j] < Jmin ) J[i][j] = Jmin;
				if( J[i][j] > Jmax ) J[i][j] = Jmax;
			}
		}
		for(int j = NE; j < N; j++){
			J[i].push_back(0.0);
			if( dice() < cEI ) J[i][j] -= JEI;
		}
	}
	for(int i = NE; i < N; i++){
		J.push_back(dvec);
		for(int j = 0; j < NE; j++){
			J[i].push_back(0.0);
			if( dice() < cIE ) J[i][j] += JIE;
		}
		for(int j = NE; j < N; j++){
			J[i].push_back(0.0);
			if( i != j && dice() < cII ) J[i][j] -= JII;
		}
	}

	for(int q = 0; q < p; q++){
		for(int i = q*NEa; i < (q+1)*NEa; i++){
			for(int j = q*NEa; j < (q+1)*NEa; j++){
				J[i][j] = 0.0;
				if( i != j && dice() < cEE ){
					J[i][j] += Jp*(1.0 + sigJ*ngn());
					if( J[i][j] < Jmin ) J[i][j] = Jmin;
					if( J[i][j] > Jmax ) J[i][j] = Jmax;
				}
			}
		}
	}
		
	return J;	
}

void simul(double alpha, double usd, double JEI, double a, double Jp, int p, int k){	
	double Jb = 0.16;
	vector< vector<double> > J = calc_J(Jp, Jb, JEI, a, p);
	calc(J,alpha,usd,JEI,a,Jp,p,k);
}

int main(int argc, char **argv){
	cout << SNE << " " << SNI << endl;
	double alpha = 0.0; //weight dependence of LTD (alpha=50.0)
	double usd = 0.0; //release probability (used = 0.05-0.5)
	double JEI = 0.0; //strngth of I-to-E connections (JEI=0.15) 
	double a = 0.0; //sparseness of embedded assemblies (a=1/p)
	double Jp = 0.0; //the mean synaptic weight within an assembly (Jp=0.3-0.7)
	int p = 0; //the number of embedded patterns (p=3,5,32)
	int k = 0; //simulation id
	if(argc > 1) alpha = atof(argv[1]);
	if(argc > 2) usd = atof(argv[2]);
	if(argc > 3) JEI = atof(argv[3]);
	if(argc > 4) a = atof(argv[4]);
	if(argc > 5) Jp = atof(argv[5]);
	if(argc > 6) p = atoi(argv[6]);
	if(argc > 7) k = atoi(argv[7]);
	cout << alpha <<" "<< usd <<" "<< JEI <<" "<< a <<" "<< Jp <<" "<< p <<" "<< k << endl;
	srand((unsigned int)(time(NULL)+k));
	simul(alpha,usd,JEI,a,Jp,p,k);
	return 0;
}
