Aim: Calculating average carrying capacity after a new mutant emerges (short-term behavior)

Main Codes:
	(1) meanKGames.cpp --- using for anayzing each game type
		input: initial_payoff number_of_samples theta  
		output: one file "data/meanK_game_..." - each line indicates each game type, information for extinction, population composition, and survival's payoff
		fprintf(fp, "%d %d %d %d %g\n", sys.game, sys.ext, sys.n-1, N, sys.P[0][0]);
	
	(2) meanKfromMutant.cpp --- calculating average carrying capacity for all game types at a given theta
		input: initial_payoff number_of_samples theta  
		output: one file "data/meanK_..." - each line indicates the average carrying capacity at a given theta (total sample average and surviving sample average)
		fprintf(fp, "%g %g %g %g %g\n", theta, mean, sqrt(var)/sqrt(Nen), sur_mean, sqrt(sur_var)/sqrt(Nsur) );
