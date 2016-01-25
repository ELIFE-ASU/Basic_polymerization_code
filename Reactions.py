### Recycling Chemistry Evolution Code Reactions File - for questions email:sara.i.walker@asu.edu

from Polymers import *
from Parameters import *

import math, random

from operator import itemgetter

############################################################################################################################################

def check_new(seq):
    """Checks new sequence against dictionary of sequences, adds +1 to dictionary if found, else adds new index to dictionary if new"""

    import Parameters

    new_seq_found = False


    for ID1, v in sequences.items():

    	"""Iterates over all objects in the sequence list"""

    	if sequences[ID1].seq == seq:

    		"""Checks new sequence against sequences in current dictionary"""

    		#update tot num. sequence if found

    		sequences[ID1].tot += 1
    		new_seq_found = True

    		break

    if new_seq_found == False:
		"""If new sequence was not found in dictionary, add it to library"""

		new_ID = Parameters.tot_species

		sequences[new_ID] = polymer(1, len(seq), seq, new_ID)

		Parameters.tot_species += 1

		
    #return new_ID

############################################################################################################################################

def polymerization(sequences):
	"""Calculate polymerization event"""

	#print 'polymerizing ...'

	import Parameters

	from Parameters import kp

	dice_roll = random.random()*Parameters.Ap_p #random draw of polymerization rxn to execute

	checkpoint = 0.0

	rxn_found = False

	for ID1, v in sequences.items():

		"""Iterates over all objects in the sequence list"""

		if rxn_found == True:

			break

		if sequences[ID1].tot > 0:

			for mID in range(N_mono_species):

				"""Calculate propensities for attachment of a monomer to growing sequence"""

				if sequences[mID].tot > 0 and ID1 != mID:
					"""Calculates propensity for attachment of monomers to growing polymer with Length >= 2, or to monomer of different species"""
					checkpoint += kp*sequences[ID1].tot*sequences[mID].tot

					#print 'check1'

				elif sequences[ID1].tot > 1 and ID1 == mID:
					"""Avoids double counting if mID == ID1"""
						
					checkpoint += kp*sequences[ID1].tot*(sequences[ID1].tot - 1)

					#print 'check2'


				if checkpoint > dice_roll:

					"""Rxn selected"""

					rxn_found = True

					new_sequence = sequences[ID1].seq + sequences[mID].seq  #attaches monomer to right end of growing string (unidirectional growth)

					#print sequences[mID].tot

					#update sequence counts for sequences lost in polyerization process
					sequences[mID].tot -= 1
					sequences[ID1].tot -= 1

					#print sequences[mID].tot

					if sequences[mID].tot < 0:

						print mID, ID1

						print 'mID has gone negative'

						exit()

					if sequences[ID1].tot < 0:

						print mID, ID1

						print 'ID1 has gone negative'

						exit()

					if sequences[ID1].length == 1:
						"""Updates Polymer and Monomer counts if both reacting sequences are monomers"""

						#Update count of total # of monomers
						Parameters.Nmono -= 2
						Parameters.Npoly += 1

					else:
						"""Update monomer count if only one reacting seq is a monomer (polymer tot count does not change)"""

						Parameters.Nmono -= 1

					break


	"""Check to see if new sequence is in current dictionary"""

	check_new(new_sequence)

	Parameters.kp_events += 1


############################################################################################################################################

def degradation(sequences):
    """Calculate degradation event"""

    #print 'degrading polymer ...'

    import Parameters

    from Parameters import kh

    dice_roll = random.random()*Parameters.Ap_h

    checkpoint = 0.0

    for ID, v in sequences.items():
    	"""Iterates over all objects in the sequence list"""

    	if sequences[ID].tot > 0 and sequences[ID].length > 1:
    		"""Calculate degradation for polymers with length L >=2"""

    		checkpoint += kh*sequences[ID].tot*(sequences[ID].length - 1) #equal probability of breaking at any bond#

    	if checkpoint > dice_roll:

    		sequences[ID].tot -= 1

    		Parameters.Npoly -= 1 #remove polymer from total count 

    		break

    break_bond(sequences, ID)

    Parameters.kh_events += 1

############################################################################################################################################

def break_bond(sequences, ID):
    """chooses bond at random in sequence for degration event"""

    import Parameters

    #print sequences[ID].seq

    b = random.randint(1, sequences[ID].length -1) #draws random number between 0 and number of bonds

    #print b
    
    #breaks sequence at bond b
    seq1 = sequences[ID].seq[:b]
    seq2 = sequences[ID].seq[b:]

    """Check if new sequences are in extant sequence population"""

    check_new(seq1)
    check_new(seq2)

    """Updates total polymer and monomer counts"""

    if len(seq1) > 1:

    	Parameters.Npoly += 1

    else:

    	Parameters.Nmono += 1

    if len(seq2) > 1:

    	Parameters.Npoly += 1

    else:

    	Parameters.Nmono += 1

######################################################################################################################
def replication(sequences):
	"""Calculates replication events (no mutation)"""

	#print 'replicating ...'

	import Parameters

	from Parameters import N_mono_species, kr

	dice_roll = random.random()*Parameters.Ap_r

	checkpoint = 0.0


	for ID, v in sequences.items():

		if sequences[ID].tot > 0 and sequences[ID].length >= R_L:

			"""Check sufficient resources for replication"""
			enough_food = True

			for mID in range(N_mono_species):

				#print 'checking mono availability'
				#print ID, mID
				#print sequences[ID].seq.count(monomer_species[mID]), sequences[mID].tot
				#print sequences[ID].seq, sequences[ID].length

				if sequences[ID].seq.count(monomer_species[mID]) > sequences[mID].tot:
					"""Compares # monomers in sequence to that available in environment"""

					enough_food = False

					break


			if enough_food == True:
				"""Calculate propensity for replication, if enough food to replicate"""

				nucleate = [int(i) for i in list(sequences[ID].seq)]

				for n in range(len(nucleate) - 1):

					mID1 = nucleate[n]
					mID2 = nucleate[n+1]

					checkpoint += kr*sequences[ID].tot*sequences[mID1].tot*sequences[mID2].tot


		if checkpoint > dice_roll:

			#update polymer count
			sequences[ID].tot += 1

			#update monomer counts

			for mID in range(N_mono_species):

				sequences[mID].tot -= sequences[ID].seq.count(monomer_species[mID])


			#update polymer and monomer counts
			Parameters.Npoly += 1
			Parameters.Nmono -= sequences[ID].length

			break


	Parameters.kr_events += 1










 


