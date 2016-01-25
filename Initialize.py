### Recycling Chemistry Evolution Code Initialization File - for questions email:sara.i.walker@asu.edu

import Polymers

###############################################################################
def perm_unique(elements):
    """Generates all unique permutations of list in elements"""
    eset=set(elements)
    listunique = [[i,elements.count(i)] for i in eset]
    u=len(elements)
    return perm_unique_helper(listunique,[0]*u,u-1)

def perm_unique_helper(cl,lst,d):
    if d < 0:
        yield tuple(lst)
    else:
        for (ind,(j,i)) in enumerate(cl):
            if i>0:
                lst[d]=j
                cl[ind][1]-=1
                for g in  perm_unique_helper(cl,lst,d-1):
                    yield g
                cl[ind][1]+=1


###############################################################################

def initialize_sequences(exp, sequences):
    """This initializes the polymer list, which usually consists only of monomers to start"""

    import Parameters

    from Parameters import monomer_species, init_monomer_count, N_mono_species

    for mID in range(0, N_mono_species):
        """Loop over monomer species in system"""

        sequences[mID] = Polymers.polymer(init_monomer_count[mID], 1, monomer_species[mID], mID)  #adds monomer to sequence list
        
        Parameters.Nmono += sequences[mID].tot #adds tot # monomers to Nmono

        Parameters.tot_species += 1 #updates total number of species count


#####################################################################################################

def initialize_Parameters():
    """Reinitializes Parameters for each new experimental run"""

    import Parameters

    Npoly = 0   		# Total number of polymers in system (summed over lattice sites)
    Nmono = 0               # Total number of monomers in system
    tot_species = 0         # Total number of replicator species that have ever existed in sim

    khevents = 0            # Tracks number of polymer replication events
    khevents = 0            # Tracks number of polymer degradation events
    kpevents = 0            # Tracks number of polymer polymerization events
   


#####################################################################################################
###########    Initialize and Update Reaction Propensities    #######################################
#####################################################################################################
def update_propensities(sequences):
	"""Initialize or update reaction propensities"""

	#print 'updating propensities ... '

	import Parameters
	from Parameters import monomer_species, N_mono_species, kp, kh, kr, R_L

	#Initialize local variables to store propensities
	Atot = 0.0
	Ap_r = 0.0
	Ap_h = 0.0
	Ap_p = 0.0

    
	for ID1, v in sequences.items():
		"""Iterates over all objects in the sequence list"""

		if sequences[ID1].tot > 0:

			for mID in range(N_mono_species):
				"""Calculate propensities for polymerization"""


				if sequences[mID].tot > 0 and ID1 != mID:

					"""Calculates propensity for attachment of monomers to growing polymer with Length >= 2, or to monomer of different species"""

					Ap_p += kp*sequences[ID1].tot*sequences[mID].tot

				elif sequences[ID1].tot > 1 and ID1 == mID:
					"""Avoids double counting if mID == ID1"""

					Ap_p += kp*sequences[ID1].tot*(sequences[ID1].tot - 1)


			if sequences[ID1].length > 1:

				"""Calculate degradation for polymers with length L >=2"""

				Ap_h += kh*sequences[ID1].tot*(sequences[ID1].length - 1) #equal probability of breaking at any bond



			if sequences[ID1].length >= R_L and kr != 0:
				"""Calculate replication for polymers with length >= R_L"""

				"""Check sufficient resources for replication"""
				enough_food = True

				for mID in range(N_mono_species):

					if sequences[ID1].seq.count(monomer_species[mID]) > sequences[mID].tot:
						"""Compares # monomers in sequence to that available in environment"""

						enough_food = False

						break


				if enough_food == True:
					"""Calculate propensity for replication, if enough food to replicate"""

					nucleate = [int(i) for i in list(sequences[ID1].seq)]

					for n in range(len(nucleate) - 1):

						mID1 = nucleate[n]
						mID2 = nucleate[n+1]

						Ap_r += kr*sequences[ID1].tot*sequences[mID1].tot*sequences[mID2].tot




	if Ap_p < 0 or Ap_h < 0 or Ap_r < 0:

		print 'Negative propensity encountered, existing sim ...'

		print Ap_p, Ap_h, Ap_r

		exit()

	#print Ap_p, Ap_h, Ap_r

	Parameters.Ap_p = Ap_p  #update propensities globally for other functions to access
	Parameters.Ap_h = Ap_h

	Parameters.Ap_r = Ap_r
	Parameters.Atot = Ap_p + Ap_h + Ap_r






   
