### Recycling Chemistry Evolution Code Main File - for questions email:sara.i.walker@asu.edu

#import pp
import sys
from time import clock
from Polymers import *
from Initialize import *
from Reactions import *
from Output import *

import numpy as np
import random
import math

def main():

    import Parameters

    for exp in range(0, Total_runs):

        start_time = clock()

        """Initialize random number generator for each run"""
        Parameters.run_seed = 100*exp + run_seed  
        random.seed(Parameters.run_seed)

        """Clear variables for each new experimental run"""
        sequences.clear()
        initialize_Parameters()
        initialize_sequences(exp, sequences)

        update_propensities(sequences)

        tau = 0 #initialize time = 0 for new run
        freq_counter = 0#counter used to keep track of frequent data outputs
        infreq_counter = 0  #counter used to keep track of infrequent data outputs
        infreq_counter = output_data(exp, infreq_counter, 0, sequences) #outputs t=0 statistics

        while tau < tau_max:
            """Main time evolution loop"""

            """Choose Reaction and Execute in Gillespie Algorithm"""

           
            dice_roll = random.random()*Parameters.Atot

            if(dice_roll < Parameters.Ap_p):
                polymerization(sequences)

            elif(dice_roll < (Parameters.Ap_p + Parameters.Ap_h)):
                degradation(sequences)

            elif(dice_roll < (Parameters.Ap_p + Parameters.Ap_h + Parameters.Ap_r)):
                replication(sequences)


            """Update rxn propensities after rxn is executed"""
            update_propensities(sequences)

            #print 'clock time ' + str(tau)

            #############################################################
            """Check Mass Conservation - If violated evolution loop will exit"""

            mass = 0

            for ID, v in sequences.items():

                if sequences[ID].tot < 0:

                    print 'seq ' + str(ID) + ' has gone negative'

                    exit()

                mass += sequences[ID].length*sequences[ID].tot

            if mass != sum(init_monomer_count):

                print 'Conservation of mass violated, exiting simulation ... '

                exit()


            #############################################################

            """Update Simulation Time (2nd part of Gillespie loop)"""

            dice_roll = random.random()
            tau -= math.log(dice_roll)/Parameters.Atot

            if tau > freq_counter:

                """Output data to file"""

                print tau

                infreq_counter = output_data(exp, infreq_counter, tau, sequences)

                freq_counter += t_frequent

        print 'Saving final data for run number ' + str(exp) + ' ...'
        runtime = clock() - start_time
        final_data(exp, runtime)
        output_data(exp, infreq_counter, tau, sequences)

# End Function Main
if __name__=="__main__":
    main()

