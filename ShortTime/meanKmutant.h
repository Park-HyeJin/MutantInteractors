#ifndef mutant_h
#define mutant_h
#include"twist.h"
#include <iostream>
#include <fenv.h>
#include <limits>
using namespace std;

struct Pop
{
	//for the system
	int game; //0: Dominance of mutant, 1: Coexistence, 2: Coordination, 3: Dominance of residence
	int n; //# of types
	int* x; //abundance vector
	double** P; //P[i][j] means the payoff of i type with an j type opponent
	double** rates; //1./P[i][j]/M
	double mu; //mutation probability
	int M; //typical population size

	//for birth and death rates
	double rd;
	double rb;
	double r; //r=rb-rd

	//for dynamics
	double tot;
	int choice;
	double theta;
	double t;
	double dt;

	//for writting
	int ns;//maximum # of types at the same time
	int id;//total number of types including extincted types
	int ext;
};


void init(Pop &sys, int M, double mu, double theta)
{
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

	sys.t =0;
	sys.dt = 0.5;
	sys.rd=0.4;
	sys.rb=0.9;
	sys.r = sys.rb-sys.rd;
}

double exponential(Pop &sys, double d)
{
	//determine lambda satisfying desired theta value
	double lambda = -log(sys.theta)/d;

	//sampling from the distribution lambda*exp(-lambda x)
	return -log(1.-drnd())/lambda;
}

void exponential_sampling(Pop &sys, int mom)
{
	for(int i=0; i<sys.n; i++)
	{
		sys.P[sys.n][i] = exponential(sys, sys.P[mom][i]);
		sys.P[i][sys.n] = exponential(sys, sys.P[i][mom]);
	}
	sys.P[sys.n][sys.n] = exponential(sys, sys.P[mom][mom]);
}

void rescale(Pop &sys)
{
	//P dividied by max and calculate reaction rates
	for(int i=0; i<sys.n; i++)
	{
		for(int j=0; j<sys.n; j++)
		{
			sys.rates[i][j] = 1./sys.P[i][j]/sys.M;
		}
	}
}

//0: Dominance of mutant, 1: Coexistence, 2: Coordination, 3: Dominance of residence
int gametype(Pop &sys)
{
	int game;
	if(sys.P[0][0] < sys.P[1][0])
	{
		if(sys.P[1][1] < sys.P[0][1])
		{
			game = 1;
		}
		else
		{
			game = 0;
		}
	}
	else
	{
		if(sys.P[1][1] < sys.P[0][1])
		{
			game = 3;
		}
		else
		{
			game = 2;
		}
	}
	return game;
}

void mutate(Pop &sys, int mom)
{
	//pick the new payoffs
	exponential_sampling(sys, mom);
	sys.x[sys.n] = 1;
	sys.n++;
	sys.id++;

	//rescale the payoffs and and calculate the rates
	rescale(sys);
}


void init_conf(Pop &sys, double a)
{
	sys.n = 1;
	sys.id = 1;
	sys.t = 0;
	sys.ext=-1;
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
	sys.P[0][0] = a;
	sys.rates[0][0] = a/(double)M;

	//mutation
	mutate(sys, 0);
	sys.game = gametype(sys);
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

	if(choice==-1){cout<< "something wring in choice: tot=" << tot << endl;  exit(1);}
	free(TNAS);
}


void extinction(Pop &sys, int die)
{
	sys.ext = die;
	sys.n -= 1;
	sys.x[die] = sys.x[sys.n];
	for(int i=0; i<sys.n; i++)
	{
		if(i==die)
		{
			sys.P[die][i] = sys.P[sys.n][sys.n];
			sys.rates[die][i] = sys.rates[sys.n][sys.n];
		}
		else
		{
			sys.P[die][i] = sys.P[sys.n][i];
			sys.P[i][die] = sys.P[i][sys.n];
			sys.rates[die][i] = sys.rates[sys.n][i];
			sys.rates[i][die] = sys.rates[i][sys.n];
		}
	}
}

void update(Pop &sys) {
	int i, j;

	//birth: choice \in [0:n-1]
	if(sys.choice < sys.n) {
		if( drnd() < sys.dt*sys.rb ) {
			//show(sys);
			if(drnd() < sys.mu)	{
				mutate(sys, sys.choice);
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
				extinction(sys, sys.choice);
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
				extinction(sys, i);
			}
		}
	}

	//update t
	sys.t += sys.dt/sys.tot;
}



void run(Pop &sys)
{
	choose_n(sys);
	update(sys);
}


#endif
