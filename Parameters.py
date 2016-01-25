### Recycling Chemistry Evolution Code Parameter File - for questions email:sara.i.walker@asu.edu

#Import Python Modules
from math import sqrt, pi, pow
import os

#####################################################################################################
run_seed = 3500
Total_runs = 10              #Number of experimental runs - used for cluster computing.


"""System Parameters & Output Specifications"""
tau_max = 500                                #max simulation time (in dimensionless time units)
t_frequent = 1.0
t_infrequent = 1                #total number of Gillespie steps between infrequent outputs
#N_bins = 10                                    #Total number of bins for population k_value histogram

print_all = 0              #If = 0 no polymer concentrations output, if = 1 polymer concentrations output for all polymers at each time step
print_plots= 1             #IF = 1, outputs plots for infrequent calculations, if 0 outputs data files


"""Replicator statistics"""
R_L = 5                     # Minimal length for replicators

"""Microscopic Reaction Rates"""
kp = 0.00001 #0.0005          #polymerization rate constant
kh = 1.0                #degradation rate constant
kr = 1.0             #replication rate constant

"""Monomer Species & Statistics"""
monomer_species = ['0', '1']

init_monomer_count = [10000, 10000]                 #Initial number of each monomeric species at start of simulation
N_mono_species = len(monomer_species)


####################################################################################################
### Do not modify parameters below this line!
#####################################################################################################

"""Global variables not to be changed (program will modify)"""

sequences = {}          #List of all sequences present in the system
monomers =  {}          #List of monomer species and abundances


Npoly = 0   		# Total number of polymers in system (summed over lattice sites)
Nmono = 0               # Total number of monomers in system
tot_species = 0         # Total number of replicator species that have ever existed in sim


Atot = 0.0                #Total global propensity for all events
Ap_r = 0.0                #Total propensity for replication
Ap_h = 0.0                #Total propensity for degradation
Ap_p = 0.0                #Total propensity for polymerization

kr_events = 0            # Tracks number of replication events
kh_events = 0            # Tracks number of degredation events
kp_events = 0            # Tracks number of polymerization events

mass = sum(init_monomer_count)

dirname = ('data/kp%.5f_kr%.5f_kh%.3f_N%d' % (kp, kr, kh, mass))

if not os.path.exists(dirname):
	"""Creat directory with dirname, if it does not exist"""
	os.makedirs(dirname)
