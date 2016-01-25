### Recycling Chemistry Evolution Code Output File - for questions email:sara.i.walker@asu.edu


import Parameters
import numpy as np
import matplotlib.pyplot as plt
#from pylab import *
import math

############################################################

def output_data(exp, infreq, tau, sequences):
    """Output run data"""

    from Parameters import t_infrequent

    """Frequent Calculations"""

    Abundances(exp, tau, sequences) #prints Npoly and Nmono to file

    output_monomers(sequences, exp, tau) #prints monomer abundances to file

    Avg_length_and_extant_species(sequences, exp, tau) #prints avg. length and total number of unique sequences to file

    
    if Parameters.print_all == 1:
    	"""Outputs all sequences in extant pop if == 1"""

        output_all(exp, tau, sequences)
        

    """Infrequent Calculations"""
    if tau > 0 and tau > infreq:


        length_distribution(sequences, exp, infreq, tau)  #length distribution of polymers

        
        infreq += t_infrequent


    return infreq


##########################################################################################################################
# Frequent Calculations
#################################################################################################################

def Abundances(exp, t, sequences):

	"""Print total number of extant polymer species as a function of time"""
	import Parameters

	"""
	poly_count = 0.0
	mono_count = 0.0

	for ID, v in sequences.items():

		if sequences[ID].length != 1:

			poly_count += sequences[ID].tot

		else:

			mono_count += sequences[ID].tot
	"""
  
	filename = ('%s/%i_tot_monomers.dat' % (Parameters.dirname, exp))

	if(t == 0):
		file = open(filename, 'w')
	else:
		file = open(filename, 'a')

	s = str(t) + '       ' + str(Parameters.Nmono)
    
	file.write(s)
	file.write('\n')
	file.close()

	filename = ('%s/%i_tot_polymers.dat' % (Parameters.dirname, exp))

	if(t == 0):
		file = open(filename, 'w')
	else:
		file = open(filename, 'a')

    
	s = str(t) + '       ' + str(Parameters.Npoly)
    
	file.write(s)
	file.write('\n')
	file.close()

######################################################################

def output_monomers(sequences, exp, t):
    """Print time-dependent polymer concentrations"""

    from Parameters import N_mono_species

    
    for mID in range(N_mono_species):

        filename = ('%s/%i_monomer_%s.dat' % (Parameters.dirname, exp, sequences[mID].seq))
    
        if(t == 0):
            file = open(filename, 'w')
        else:
            file = open(filename, 'a')
            
        s = str(t) + '      ' + str(sequences[mID].tot)
        file.write(s)
        file.write('\n')
        file.close()


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

    file_name = ('%s/%i_avg_L.dat' % (Parameters.dirname, exp))

    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
        
    s = str(t) + "        " + str(avg_L)
    file.write(s)
    file.write('\n')

    file.close()   

    """Write Extant Species Count to File"""

    file_name = ('%s/%i_extant_species.dat' % (Parameters.dirname, exp))
    
    if(t == 0):
        file = open(file_name, 'w')
    else:
        file = open(file_name, 'a')
            
    s = str(t) + '      ' + str(num_extant)
    file.write(s)
    file.write('\n')
    file.close() 


    
##########################################################################################################################
# Infrequent Calculations
#################################################################################################################

def length_distribution(sequences, exp, t, tau):

	#bin polymers by length up to length max_L, all seq length > max_L binned in max_L
	max_L = 10

	bin_x = max_L*[0]  
	bin_y = max_L*[0]    

	for (ID, seq) in sequences.items():

		L = seq.length - 1

		if L > 9:

			L = 9

		bin_x[L] = seq.length
		bin_y[L] += seq.tot

	if Parameters.print_plots == 1:

		fig = plt.figure()   
		ax = fig.add_subplot(111)
		ax.set_ylabel('Abundance')
		ax.set_xlabel('Length')
		ax.set_xticks(bin_x)
		width = 0.25

		rects1 = ax.bar(bin_x, bin_y , width, color='blue', align = 'center')

		file_name = ('%s/%i_Length_Dist_%i.png' % (Parameters.dirname, exp, t))
		plt.savefig(file_name, format = 'png')

	elif Parameters.print_plots == 0:

		#Save length distribution to a numpy array   
		file_name = ('%s/%i_Length_Dist_%i' % (Parameters.dirname, exp, t))
		np.save(file_name, bin_y)



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

	file.write("Number of replication events = ")
	s = str(Parameters.kr_events)
	file.write(s)
	file.write('\n')

	file.write("Number of degredation events = ")
	s = str(Parameters.kh_events)
	file.write(s)
	file.write('\n')

	file.close()


################################################################################



