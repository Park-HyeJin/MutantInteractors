#ifndef mutant_h
#define mutant_h
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include"twist.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <fenv.h>
using namespace std;


struct Pop {
	//for the system
	int n; //# of types
	int* x; //abundance vector
	double** P; //P[i][j] means the payoff of i type with an j type opponent
	double** rates; //1./P[i][j]/M
	double mu; //mutation probability
	int M; //typical population size
	int ns;//maximum # of types at the same time
	int id;//total number of types including extincted types; id of the first residence is id=0

	//for birth and death rates
	double rd; //death rate
	double rb; //birth rate
	double r; //r=rb-rd

	//for dynamics
	double tot; //Total number of accessible states
	int choice; //reaction idx which will occur
	double theta; //prob for advantageous mutation
	double t; //real time
	double dt; //t += dt/tot
	int maxt;//unit of mutant (generations); is the same with maximum number of id
};

struct Data {
	//for writting
	int maxt;
	int Nt; //the number of recording time
	double delt;
	int fidx;
	double* tempt;
	double* tempn;
	double* tempN;
	double** tempx;
	double** Payoffs;
};

void extinction(Pop &sys, Data &mut, int die) {
	sys.n -= 1;
	sys.x[die] = sys.x[sys.n];
	for(int i=0; i<sys.n; i++) {
		if(i==die) {
			sys.P[die][i] = sys.P[sys.n][sys.n];
			sys.rates[die][i] = sys.rates[sys.n][sys.n];
		}
		else {
			sys.P[die][i] = sys.P[sys.n][i];
			sys.P[i][die] = sys.P[i][sys.n];
			sys.rates[die][i] = sys.rates[sys.n][i];
			sys.rates[i][die] = sys.rates[i][sys.n];
		}
	}

	//record extinction of whole population event
	if(sys.n==0) {
		mut.tempn[sys.id] = 0;
		mut.tempN[sys.id] = 0;
		mut.Payoffs[sys.id++][0] = sys.P[0][0];
	}
}

double exponential(Pop &sys, double d) {
	//determine lambda satisfying desired theta value
	double lambda = -log(sys.theta)/d;

	//sampling from the distribution lambda*exp(-lambda x)
	double res = -log(1.-drnd())/lambda;
//	cout << res <<  " ";
	return res;
}

void exponential_sampling(Pop &sys, int mom) {
	for(int i=0; i<sys.n; i++) {
		sys.P[sys.n][i] = exponential(sys, sys.P[mom][i]);
		sys.P[i][sys.n] = exponential(sys, sys.P[i][mom]);
	}
	sys.P[sys.n][sys.n] = exponential(sys, sys.P[mom][mom]);
}

void rates(Pop &sys) {
	for(int i=0; i<sys.n; i++) 	{
		for(int j=0; j<sys.n; j++) 	{
			sys.rates[i][j] = 1./sys.P[i][j]/sys.M;
		}
	}
}

void record_mut(Pop &sys, Data &mut){
	mut.tempt[sys.id] = sys.t;
	mut.tempn[sys.id] = sys.n;
	mut.tempN[sys.id] = 0;
	for(int i=0; i<sys.n; i++) {
		mut.tempN[sys.id] += sys.x[i];
		mut.tempx[sys.id][i] =  sys.x[i];
		for(int j=0; j<sys.n; j++) { 
			mut.Payoffs[sys.id][sys.n*i+j] = sys.P[i][j];
		}
	}
}

void record_new_payoff(Pop &sys, Data &mut){
	int n0 = sys.n*sys.n;
	int n1 = sys.n+1;
	for(int i=0; i<n1; i++) {
		for(int j=0; j<n1; j++) { 
			mut.Payoffs[sys.id-1][n0+n1*i+j] = sys.P[i][j];
		}
	}
}

void mutate(Pop &sys, Data &mut, int mom) {
	record_mut(sys, mut);
	sys.id++;

	//pick the new payoffs
	exponential_sampling(sys, mom);
	record_new_payoff(sys, mut);
	sys.x[sys.n] = 1;
	sys.n++;

	//rescale the payoffs and and calculate the rates
	rates(sys);
}

void choose_n(Pop &sys) {
	int N = (sys.n+2)*sys.n;
	double* TNAS = (double*)calloc(N, sizeof(double));
	int idx=0;
	double tot=0;

	//make prob array
	//for birth
	for(int i=0; i<sys.n; i++) {
		tot += sys.x[i];
		TNAS[idx++] = tot;
	}

	//for death
	for(int i=0; i<sys.n; i++) {
		tot += sys.x[i];
		TNAS[idx++] = tot;
	}

	//interaction between different type
	for(int i=0; i<sys.n; i++) {
		for(int j=0; j<sys.n; j++) {
			if(i==j) {
				//interaction with the same type
				tot += (double)sys.x[i]*(sys.x[i]-1);
				TNAS[idx++] = tot;
			}
			else {
				tot += (double)sys.x[i]*sys.x[j];
				TNAS[idx++] = tot;
			}
		}
	}

	//find reaction rule
	int choice=-1;
	double target = drnd()*tot;
	if(target < TNAS[0]) 	{
		choice = 0;
	}
	else 	{
		for(int i=1; i<idx; i++) {
			if(target < TNAS[i]) {
				choice = i;
				break;
			}
		}
	}
	sys.tot = tot;
	sys.choice = choice;

	if(choice==-1){cout<< "something wring in choice: tot=" << tot << endl; exit(1);}
	free(TNAS);
}

void update(Pop &sys, Data &mut) {
	int i, j;

	//birth: choice \in [0:n-1]
	if(sys.choice < sys.n) {
		if( drnd() < sys.dt*sys.rb ) {
			//show(sys);
			if(drnd() < sys.mu)	{
				mutate(sys, mut, sys.choice);
			}
			else {
				sys.x[sys.choice]++;
			}
		}
	}
	//death: choice \in [n:2n-1]
	else if(sys.choice < 2*sys.n) {
		if (drnd() <sys.dt*sys.rd )	{
			//show(sys);
			sys.choice -= sys.n; //choice \in [0:n*n-1]
			sys.x[sys.choice] -= 1;
			if(sys.x[sys.choice] == 0) {
				extinction(sys, mut, sys.choice);
			}
		}
	}
	//death from competition: choice \in [2n:n*n+2n]
	else {
		//show(sys);
		sys.choice -= 2*sys.n; //choice \in [0:n*n-1]
		i = sys.choice%sys.n;
		j = (sys.choice - i)/sys.n;
		//cout << sys.rates[i][j] << endl;

		if( drnd() < sys.dt*sys.rates[i][j] ) {
			sys.x[i] -= 1;
			if(sys.x[i] == 0) {
				extinction(sys, mut, i);
			}
		}
	}

	//update t
	sys.t += sys.dt/sys.tot;
}

void run(Pop &sys, Data &mut) {
	choose_n(sys);
	update(sys, mut);
}

void init(Pop &sys, Data &mut, int M, double mu, double theta, int maxt, int fidx) {
	//for Pop
	sys.M = M;
	sys.ns = 50;
	sys.mu = mu;
	sys.theta = theta;
	sys.x = (int*)calloc(sys.ns, sizeof(int));
	sys.P = (double**)calloc(sys.ns, sizeof(double*));
	sys.rates = (double**)calloc(sys.ns, sizeof(double**));
	for(int i=0; i<sys.ns; i++)
	{
		sys.P[i] = (double*)calloc(sys.ns, sizeof(double));
		sys.rates[i] = (double*)calloc(sys.ns, sizeof(double));
	}
	sys.maxt = maxt;
	sys.t =0;
	sys.dt=0.5;
	sys.rd=0.4;
	sys.rb=0.9;
	sys.r = sys.rb-sys.rd;

	//for mut Data structure
	mut.fidx = fidx;
	mut.maxt = maxt;
	mut.delt = 1;
	mut.Nt = maxt;
	mut.tempt = (double*)calloc(mut.Nt, sizeof(double));
	mut.tempn = (double*)calloc(mut.Nt, sizeof(double));
	mut.tempN = (double*)calloc(mut.Nt, sizeof(double));
	mut.tempx = (double**)calloc(mut.Nt, sizeof(double*));
	mut.Payoffs = (double**)calloc(mut.Nt, sizeof(double*));
	for(int i=0; i<mut.Nt; i++)
	{
		mut.tempx[i] = (double*)calloc(sys.ns, sizeof(double));
		mut.Payoffs[i] = (double*)calloc(sys.ns*sys.ns, sizeof(double));
	}
}

void init_conf(Pop &sys, Data &mut) {
	sys.n = 1;
	sys.id = 0;
	sys.t = 0;
	int M = sys.M;

	for(int i=0; i<sys.ns; i++)
	{
		sys.x[i] = 0;
		for(int j=0; j<sys.ns; j++)
		{
			sys.P[i][j] = 0;
			sys.rates[i][j] = 0;
		}
	}
	sys.x[0] = (int)(sys.r*M);
	sys.P[0][0] = 1.;
	sys.rates[0][0] = 1./(double)M;

	//intialze data array
	for(int i=0; i<mut.Nt; i++) {
		mut.tempt[i] = 0;
		mut.tempn[i] = 0;
		mut.tempN[i] = 0;
		for(int j=0; j<sys.ns*sys.ns; j++) {
			mut.Payoffs[i][j] = 0;
		}
		for(int j=0; j<sys.ns; j++) {
			mut.tempx[i][j] = 0;
		}
	}

	//mutation; starting from the emergence of a mutant
	mutate(sys, mut, 0);
}

void write_conf(Pop &sys, Data &mut) {
	int n;
	char fname[200];
	FILE *fp;
	sprintf(fname, "Mut_%g/M%d_mu%g_maxt%d_%d.d", sys.theta, sys.M, sys.mu, sys.maxt-1, mut.fidx);
	

	//write
	ofstream ofp(fname);
	for(int i=0; i<sys.id; i++)
	{
		ofp << mut.delt*i << "\t" << mut.tempt[i]<< "\t" << mut.tempN[i] << "\t" << mut.tempn[i] << "\t";
		
		//print x[i] and P[i][j]
		n = mut.tempn[i];if(n==0) {n = 1;}
		for(int l=0; l<n; l++) {
			ofp << mut.tempx[i][l] << "\t";
		}
		for(int l=0; l<2*n*n+2*n+1; l++) {
			ofp << mut.Payoffs[i][l] << "\t";
		}
		ofp << endl;
	}
}

#endif
