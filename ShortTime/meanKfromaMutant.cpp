/////////////////////////////////////////////////
// initial condition : one type with payoff a = 1
// output:

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "twist.h"
#include "meanKmutant.h"
using namespace std;

int main(int argc, char** argv)
{
	if(argc!=4)
	{
		printf("Usage: %s a Nen theta\n", argv[0]);
		exit(1);
	}

	double K=0;
	double K2=0;
	int N;
	int M = 1000;
	double mu = 1e-5;
	double a = atof(argv[1]);
	int Nen = atoi(argv[2]);
	int Nsur = Nen;
	double theta = atof(argv[3]);

	init_rnd(gus()+theta);
	//initialize
	Pop sys;
	init(sys, M, mu, theta);

	//run: terminate when the second mutant emerges or the population goes to extinction
	for(int en=0; en<Nen; en++)
	{
		init_conf(sys, a);
		while(sys.id==2 && sys.n!=0)
		{
			run(sys);
		}
		if(sys.n==0)//extinction
		{
			Nsur -= 1;
		}
		
		//calculating the population size
		N = 0;
		for(int i=0; i<sys.n-1; i++)
		{
			N += sys.x[i];
		}
		K += N;
		K2 += N*N;
	}

	//write the file: M r Nen K K2
	FILE *fp;
	char fname[200];
	sprintf(fname, "data/meanK_a%g_Nen%d.d", a, Nen);
	fp = fopen(fname, "a");
	double mean = K/Nen;
	double var = K2/Nen -  mean*mean;
	double sur_mean = K/Nsur;
	double sur_var = K2/Nsur -  sur_mean*sur_mean;
	fprintf(fp, "%g %g %g %g %g\n", theta, mean, sqrt(var)/sqrt(Nen), sur_mean, sqrt(sur_var)/sqrt(Nsur) );

	return 0;
}
