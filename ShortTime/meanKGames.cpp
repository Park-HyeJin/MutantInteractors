/////////////////////////////////////////////////
// initial condition : one type with payoff a=1

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

	int M = 1000;
	double mu = 1e-5;
	double a = atof(argv[1]);//initial payoff
	int Nen = atoi(argv[2]);//number of samples
	int Nsur = Nen;//number of surviving samples
	double theta = atof(argv[3]);
	int N;//population size

	FILE *fp;
	char fname[200];
	sprintf(fname, "data/meanK_game_a%g_theta%g_Nen%d.d", a, theta, Nen);
	fp = fopen(fname, "w");
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
		//recoed game type and population size
		//sys.ext=-1: both coexists
		//sys.ext=0: residents go to extinction
		//sys.ext=1: mutants go to extinction
		fprintf(fp, "%d %d %d %d %g\n", sys.game, sys.ext, sys.n-1, N, sys.P[0][0]);
	}

	fclose(fp);

	return 0;
}

