# ----------------------------------------
# USAGE:

# python res_res_distances.py config_file

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from distance_functions import *

zeros = np.zeros
square = np.square
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# VARIABLE DECLARATION

config_file = sys.argv[1]

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

necessary_parameters = ['pdb','important','traj_loc','start','end','avg_dist_output','std_dist_output']
all_parameters = ['pdb','important','traj_loc','start','end','avg_dist_output','std_dist_output','write_summary','summary_filename']
def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['write_summary'] = False 
	parameters['summary_filename'] = 'avg_structure.summary'

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

def summary():
	with open('%s' %(parameters['summary_filename']),'w') as f:
		f.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		f.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			f.write('%s ' %(sys.argv[i]))
		f.write('\n\n')
		f.write('Parameters used:\n')
		for i in all_parameters:
			f.write('%s = %s \n' %(i,parameters[i]))
		f.write('\n\n')

# ----------------------------------------
# MAIN PROGRAM:
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

start = parameters['start']
end = parameters['end']

# ----------------------------------------
# LOAD IN THE UNIVERSE TO BE ANALYZED AND INITIATE THE NECESSARY ATOM SELECTIONS
u = MDAnalysis.Universe(parameters['pdb'])
u_important = u.select_atoms(parameters['important'])

nNodes = len(u_important.residues)

# ----------------------------------------
# ARRAY DECLARATION
avg_matrix = zeros((nNodes,nNodes))
std_matrix = zeros((nNodes,nNodes))

# ----------------------------------------
# Trajectory Analysis 
nSteps = 0
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new(parameters['traj_loc'] %(start,start))
	nSteps += len(u.trajectory)

	for ts in u.trajectory:
		if ts.frame%1000 == 0:
			ffprint('Working on timestep %d of trajectory %d' %(ts.frame, start))

		for i in range(nNodes-1):
			res0 = u_important.residues[i]
			com0 = res0.center_of_mass()
			for j in range(i+1,nNodes):
				res1 = u_important.residues[j]
				com1 = res1.center_of_mass()
				dist, dist2 = euclid_dist(com0,com1)
				avg_matrix[i,j] += dist
				std_matrix[i,j] += dist2
	start +=1

# ----------------------------------------
# Finishing Averages
avg_matrix /= nSteps
std_matrix /= nSteps
std_matrix = sqrt(std_matrix - square(avg_matrix))

# ----------------------------------------
# Outputting avg and std dev arrays to files
with open(parameters['avg_dist_output'],'w') as W, open(parameters['std_dist_output'],'w') as Y:
	for i in range(nNodes):
		for j in range(nNodes):
			W.write('%10f   ' %(avg_matrix[i,j]))
			Y.write('%10f   ' %(std_matrix[i,j]))
		X.write('\n')
		Y.write('\n')

# ----------------------------------------
# Finishing Averages
if parameters['write_summary']:
	summary()

