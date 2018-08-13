Aim: Getting trajectory for each sample

* Before executing the simulation code, you needs to create a folder "Mut_theta" where theta is a parameter value which you want to use.

Main Code: singlerun.cpp
	input: M mu theta maxt(maximum mutation event time) fidx(file_index)
	output: record the trajectory of the population composition with a file index
		file_name: "Mut_theta/...fidx.d"
		each line records followings (see function "write_conf" in mutant.h file):
		mutation_event_time corresponding_real_time population_size number_of_types abundances current_payoffs new_payoffs



