#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "twist.h"
#include "mutant.h"
using namespace std;


int main(int argc, char** argv)
{
	if(argc!=6)
	{
		printf("Usage: %s M mu theta maxt(in mutantation event time) fidx\n", argv[0]);
		exit(1);
	}
	
	int M = atoi(argv[1]);
	double mu = atof(argv[2]);
	double theta = atof(argv[3]);
	int maxt = atoi(argv[4])+1;
	int fidx = atoi(argv[5]);

	//initialize
	Pop sys;
	Data mut;
	init(sys, mut,  M, mu, theta, maxt, fidx);
	init_rnd(gus()+M*theta*mu*M+fidx);
	init_conf(sys, mut);

	//run
	while(sys.id<maxt && sys.n!=0)	{
		run(sys, mut);
	}
	write_conf(sys, mut);

	return 0;
}
