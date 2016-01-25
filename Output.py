### Recycling Chemistry Evolution Code Output File - for questions email:sara.i.walker@asu.edu


import Parameters
import numpy as np
import matplotlib.pyplot as plt
#from pylab import *
import math
import csv

############################################################

def output_data(exp, infreq, tau, sequences):
	"""Output run data"""

	from Parameters import t_infrequent

	"""Frequent Calculations"""

	Abundances(exp, tau, sequences) #prints Npoly and Nmono to file

	output_monomers(sequences, exp, tau) #prints monomer abundances to file

	Avg_length_and_extant_species(sequences, exp, tau) #prints avg. length and total number of unique sequences to file


	"""Infrequent Calculations"""
	if tau > 0 and tau > infreq:

		length_distribution(sequences, exp, infreq, tau)  #length distribution of polymers and diversity measures
		
		infreq += t_infrequent


	return infreq


##########################################################################################################################
# Frequent Calculations
#################################################################################################################

def Abundances(exp, t, sequences):

	"""Print total number of extant polymer species as a function of time"""
	import Parameters

  
	filename = ('%s/%i_tot_monomers.dat' % (Parameters.dirname, exp))

	if(t == 0):
		mono_file = open(filename, 'w')
	else:
		mono_file = open(filename, 'a')

	
	mono_out = csv.writer(mono_file)
	
	mono_out.writerow((t, Parameters.Nmono))

	mono_file.close()



	filename = ('%s/%i_tot_polymers.dat' % (Parameters.dirname, exp))

	if(t == 0):
		poly_file = open(filename, 'w')
	else:
		poly_file = open(filename, 'a')


	poly_out = csv.writer(poly_file)
	
	poly_out.writerow((t, Parameters.Npoly))

	poly_file.close()


######################################################################

def output_monomers(sequences, exp, t):
	"""Print time-dependent polymer concentrations"""

	from Parameters import N_mono_species

	
	for mID in range(N_mono_species):

		filename = ('%s/%i_monomer_%s.dat' % (Parameters.dirname, exp, sequences[mID].seq))
	
		if(t == 0):
			mono_file = open(filename, 'w')
		else:
			mono_file = open(filename, 'a')
		

		mono_out = csv.writer(mono_file)
	
		mono_out.writerow((t, sequences[mID].tot))

		mono_file.close()

################################################################### 
   
def Avg_length_and_extant_species(sequences, exp, t):
	"""Outputs average length of extant sequences and num of extant species to file"""

	import Parameters

	sum_L = 0.0

	num_extant = 0


	for ID, v in sequences.items():

		sum_L += sequences[ID].length*sequences[ID].tot

		if sequences[ID].tot != 0:

			num_extant += 1

	

	avg_L = sum_L/(Parameters.Npoly + Parameters.Nmono)  #divide number of sequences by total mass to get average length

	"""Write Average Length to File"""	

	file_name = ('%s/%i_avg_L.csv' % (Parameters.dirname, exp))

	if(t == 0):
		avgL_file = open(file_name, 'w')
	else:
		avgL_file = open(file_name, 'a')


	avgL_out = csv.writer(avgL_file)
	
	avgL_out.writerow((t, avg_L))

	avgL_file.close()
	


	"""Write Extant Species Count to File"""

	file_name = ('%s/%i_extant_species.csv' % (Parameters.dirname, exp))
	
	if(t == 0):
		es_file = open(file_name, 'w')
	else:
		es_file = open(file_name, 'a')


	es_out = csv.writer(es_file)
	
	es_out.writerow((t, num_extant))

	es_file.close()
	
##########################################################################################################################
# Infrequent Calculations
#################################################################################################################

def length_distribution(sequences, exp, t, tau):

	#bin polymers by length up to length max_L, all seq length > max_L binned in max_L
	max_L = 21

	length_dist = max_L*[0]  


	GD = 0  # tracks global shannon diversity
	SD_avg = 0  # tracks average sequence diversity


	if t == 0:

		LD_file = open('%s/%i_length_dist.csv' % (Parameters.dirname, exp), 'w')

		lengths = [i for i in range(max_L)]
		
		lengths[0] = 'time'  #first index in array of lengths is the current time

		LD_out = csv.writer(LD_file)

		LD_out.writerow(lengths) 


	else:
		LD_file = open('%s/%i_length_dist.csv' % (Parameters.dirname, exp), 'a')
	

	LD_out = csv.writer(LD_file)

	length_dist[0] = tau #initilizes first bin in array to hold the time
	

	for (ID, seq) in sequences.items():

		L = seq.length

		if L >= max_L:

			L = max_L - 1  #bins all seq. of length L > 20 as L = 20

		length_dist[L] += seq.tot


		if sequences[ID].seq.count('0') != 0:

			p_i = float(sequences[ID].seq.count('0'))/float(sequences[ID].length)

        	SD_avg -= p_i*math.log(p_i, 2)



	"""Write length distribution to file"""
	LD_out.writerow(length_dist) 
	LD_file.close()
	
	Tot_pop_size = sum(length_dist[1:])

	"""Calculate and Output Global Population Diversity"""
	for (ID, seq) in sequences.items():

		if sequences[ID].tot != 0 and Tot_pop_size != 0:

			p_i = float(sequences[ID].tot)/float(Tot_pop_size)

			GD -= p_i*math.log(p_i, 2)


	if t == 0:
		GD_file = open('%s/%i_global_diversity.csv' % (Parameters.dirname, exp), 'w')

	else:
		GD_file = open('%s/%i_global_diversity.csv' % (Parameters.dirname, exp), 'a')
	


	GD_out = csv.writer(GD_file)
	GD_out.writerow((t, GD))
	GD_file.close()


	"""Calculate and Output Average sequence Diversity"""

	if Tot_pop_size != 0:
		
		SD_avg /= Tot_pop_size


	if t == 0:
		SD_file = open('%s/%i_seq_diversity.csv' % (Parameters.dirname, exp), 'w')

	else:
		SD_file = open('%s/%i_seq_diversity.csv' % (Parameters.dirname, exp), 'a')
	

	SD_out = csv.writer(SD_file)
	SD_out.writerow((t, SD_avg))
	SD_file.close()


################################################################################################################


def output_rep_abundances(exp, sequences, t):

	import Parameters
	from Parameters import init_monomer_count, T_IDs, U_IDs

	T_tot = 0.0
	U_tot = 0.0

	I_T = 0.0     #Global Diversity of ALL trivial replicator strings
	I_U = 0.0    #Global Diversity of ALL nontrivial replicator strings

	S_T_avg = 0.0     #Avg. seq diversity of all trivial replicating strings
	S_U_avg = 0.0    #Avg. seq diversity of all nontrivial replicating strings

	I_UnU = 0.0 #Diversity of Nontrivial replicators minus NT part.

	if Parameters.Trivial_Rep == True:

		for index, ID in enumerate(T_IDs):
			"""Calculates total number of trivial replicators"""

			T_tot += sequences[ID].tot


		for index, ID in enumerate(T_IDs):
			"""Calculates total diversity of trivial replicators"""

			if sequences[ID].tot != 0:

				p_i = float(sequences[ID].tot)/float(T_tot)

				I_T -= p_i*math.log(p_i, 2)


				for index, motif in enumerate(sequences[ID].T_motifs):

					parts = sequences[ID].seq.partition(sequences[motif].seq)

					nonT_seq = parts[0] + parts[2]

					if len(nonT_seq) != 0:

						p_i = float(nonT_seq.count('0'))/float(len(nonT_seq))

						if p_i != 0:

							S_T_avg -= p_i*math.log(p_i, 2)*sequences[ID].tot


		if T_tot != 0:

			S_T_avg /= T_tot  #calculates average seq. diversity of all trivial reps

		else:
			S_T_avg = 0

		T_tot /= sum(init_monomer_count)  #fraction of mass in T-replicators


		file_name = ('%s/%i_T_mass.dat' % (Parameters.dirname, exp))
	
		if(t == 0):
			Tm_file = open(file_name, 'w')
		else:
			Tm_file = open(file_name, 'a')
			

		Tm_out = csv.writer(Tm_file)
		Tm_out.writerow((t, T_tot))
		Tm_file.close()


		file_name = ('%s/%i_T_diversity.dat' % (Parameters.dirname, exp))
	
		if(t == 0):
			Td_file = open(file_name, 'w')
		else:
			Td_file = open(file_name, 'a')


		Td_out = csv.writer(Td_file)
		Td_out.writerow((t, I_T))
		Td_file.close()


		file_name = ('%s/%i_T_seq_shannon.dat' % (Parameters.dirname, exp))
	
		if(t == 0):
			Ts_file = open(file_name, 'w')
		else:
			Ts_file = open(file_name, 'a')
		

		Ts_out = csv.writer(Ts_file)
		Ts_out.writerow((t, S_T_avg))
		Ts_file.close()


	if Parameters.Universal_Rep == True:

		for index, ID in enumerate(U_IDs):
			"""Calculates total number of nontrivial replicators"""

			U_tot += sequences[ID].tot
		
		for index, ID in enumerate(U_IDs):
			"""Calculates total diversity of nontrivial replicators"""

			if sequences[ID].tot != 0:

				p_i = float(sequences[ID].tot)/float(U_tot)

				I_U -= p_i*math.log(p_i, 2)
				

				for index, motif in enumerate(sequences[ID].U_motifs):

					#print sequences[ID].seq

					parts = sequences[ID].seq.partition(sequences[motif].seq)

					#print parts
				
					#print parts[0], parts[1], parts[2]


					nonU_seq = parts[0] + parts[2]

					#print 'Non NT seq ' + str(nonNT_seq)

					if len(nonU_seq) != 0:

						p_i = float(nonU_seq.count('0'))/float(len(nonU_seq))

						if p_i != 0:

							S_U_avg -= p_i*math.log(p_i, 2)*sequences[ID].tot
						
		#exit()

		if U_tot != 0:

			S_U_avg /= U_tot  #calculates average seq. diversity of all nontrivial reps

		else:
			S_U_avg = 0

		U_tot /= sum(init_monomer_count)  #fraction of mass in T-replicators

		

		file_name = ('%s/%i_U_mass.dat' % (Parameters.dirname, exp))
	
		if(t == 0):
			Um_file = open(file_name, 'w')
		else:
			Um_file = open(file_name, 'a')


		Um_out = csv.writer(Um_file)
		Um_out.writerow((t, U_tot))
		Um_file.close()


		file_name = ('%s/%i_U_diversity.dat' % (Parameters.dirname, exp))
	
		if(t == 0):
			Ud_file = open(file_name, 'w')
		else:
			Ud_file = open(file_name, 'a')


		Ud_out = csv.writer(Ud_file)
		Ud_out.writerow((t, I_U))
		Ud_file.close()


		file_name = ('%s/%i_U_seq_shannon.dat' % (Parameters.dirname, exp))
	
		if(t == 0):
			Us_file = open(file_name, 'w')
		else:
			Us_file = open(file_name, 'a')


		Us_out = csv.writer(Us_file)
		Us_out.writerow((t, S_U_avg))
		Us_file.close()

	
################################################################################################################
		
def final_data(exp, runtime):
	"""Output final run statistics"""
	import Parameters

	"""Write Polymer list to file"""
	file = open('%s/%i_run_statistics.txt' % (Parameters.dirname, exp), 'w')

	file.write("Run parameters \n")
	s = 'Rate Constants \n'
	file.write(s)
	s = 'kp = ' + str(Parameters.kp) + '\n'
	file.write(s)
	s = 'kh = ' + str(Parameters.kh) + '\n'
	file.write(s)
	s = 'kr = ' + str(Parameters.kr) + '\n'
	file.write(s)
	 
	
	file.write("Information about kMC run including run statistics: \n \n")

	file.write("Total run time = ")
	s= str(runtime) + ' sec'
	file.write(s)
	file.write('\n')
	file.write('\n')

	file.write("Number of polymerization events")
	s = str(Parameters.kp_events)
	file.write(s)
	file.write('\n')

	file.write("Number of degradation events = ")
	s = str(Parameters.kh_events)
	file.write(s)
	file.write('\n')

	file.write("Number of replication events = ")
	s = str(Parameters.kr_events)
	file.write(s)
	file.write('\n')

	file.close()


################################################################################



