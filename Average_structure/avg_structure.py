# ----------------------------------------
# USAGE:

# python avg_structure.py config_file

# ----------------------------------------
# PREAMBLE:

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

# ----------------------------------------
# VARIABLE DECLARATION

config_file = sys.argv[1]

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

necessary_parameters = ['pdb','important','traj_loc','start','end','avg_output_pdb']
all_parameters = ['pdb','important','traj_loc','start','end','avg_output_pdb','alignment','Wrapped','substrate','write_summary','summary_filename','Convergence_Threshold','maxIterations','out_file']
def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment'] = 'protein'
        parameters['Wrapped'] = True
        parameters['substrate'] = None
	parameters['write_summary'] = False 
	parameters['summary_filename'] = 'avg_structure.summary'
        parameters['Convergence_Threshold'] = 1E-5
        parameters['maxIterations'] = 100
        parameters['out_file'] = 'average_list.txt'

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
                f.write('Average structure is written in pdb format:\n')
		f.write('	%s\n' %(parameters['avg_output)pdb']))

# ----------------------------------------
# MAIN:
# ----------------------------------------
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

start = parameters['start']
end = parameters['end']
thresh = parameters['Convergence_Threshold']
maxIter = parameters['maxIterations']

# ----------------------------------------
# LOAD IN THE UNIVERSE TO BE ANALYZED AND INITIATE THE NECESSARY ATOM SELECTIONS
u = MDAnalysis.Universe(parameters['pdb'])
u_all = u.select_atoms('all')
u_align = u.select_atoms(parameters['alignment'])
u_important = u.select_atoms(parameters['important'])
if not parameters['Wrapped']:
        u_substrate = u.select_atoms(parameters['substrate'])

# ----------------------------------------
# GRABBING IMPORTANT NUMBERS FROM THE UNIVERSE
u_align_pos = u_align.positions
u_important_atoms = u_important.n_atoms
u_align_atoms = u_align.n_atoms
if not parameters['Wrapped']:
        u_substrate_res = u_substrate.n_residues

# ----------------------------------------
# DETERMINING THE NUMBER OF STEPS TO BE AVERAGED OVER
temp = parameters['start']
nSteps = 0
while temp <= end:
	u.load_new(parameters['traj_loc'] %(temp,temp))
	nSteps += len(u.trajectory)
	temp += 1

ffprint('Number of steps to be averaged over : %d' %(nSteps))

# ----------------------------------------
# ARRAY DECLARATION
all_coord = zeros((nSteps,u_important_atoms,3),dtype=np.float32)
avgCoord = zeros((u_important_atoms,3),dtype=np.float32)
all_align = zeros((nSteps,u_align_atoms,3),dtype=np.float32)
avgAlign = zeros((u_align_atoms,3),dtype=np.float32)

# ----------------------------------------
# Trajectory Analysis 
ffprint('Beginning trajectory analysis')
temp = 0 
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new(parameters['traj_loc'] %(start,start))

	for ts in u.trajectory:
                # Removing translational motion of the system
		u_all.translate(-u_align.center_of_mass())
		
		# CALCULATIONS that are unnecessary if the trajectory is wrapped.
		if not parameters['Wrapped']:		# Test to see if the 'Wrapped' key is equal to False
			dims = u.dimensions[:3]		
			dims2 = dims/2.0
	
			for i in range(u_substrate_res):
				COM = u_substrate.residues[i].center_of_mass()
				t = wrapping(COM,dims,dims2)
				u_substrate.residues[i].atoms.translate(t)
                
#                # Removing rotational motion of the system
#                R, rmsd = rotation_matrix(u_align.positions,u_align_pos)
#                u_important.rotate(R)

                # Collecting the positions of the important and alignment selections
		avgCoord += u_important.positions
		avgAlign += u_align.positions
		all_coord[temp] = u_important.positions
		all_align[temp] = u_align.positions
		temp += 1
	start += 1

ffprint(nSteps)
if temp != nSteps:
	ffprint('Failed to collect positions from all timesteps.')
        sys.exit()

# ----------------------------------------
# Finishing Averages
avgCoord /= float(nSteps)
avgAlign /= float(nSteps)
ffprint('Finished with the trajectory analysis')

# ----------------------------------------
# Calculating and Aligning to the average positions
iteration = 0
residual = thresh + 100.0 					# arbitrary assignment greater than thresh
ffprint('Beginning iterative process of calculating average positions and aligning to the average')
while residual > thresh and iteration < maxIter:		
	tempAvgCoord = zeros((u_important_atoms,3),dtype=np.float32)		# zeroing out the tempAvgCoord array every iteration
	tempAvgAlign = zeros((u_align_atoms,3),dtype=np.float32)                # zeroing out the tempAvgAlign array every iteration
        # Looping through every step in trajectory, rotating frame to the alignment landmark, and added aligned coordinates to the temp arrays
	for i in range(nSteps):
		R, d = rotation_matrix(all_align[i,:,:],avgAlign)
#                print R
		all_align[i,:,:] = dot_prod(all_align[i,:,:],R.T)
		all_coord[i,:,:] = dot_prod(all_coord[i,:,:],R.T)
		tempAvgAlign += all_align[i,:,:]
		tempAvgCoord += all_coord[i,:,:]			# recalculate the average coordinates to optimize the average position
	# Finishing the averages
        tempAvgCoord /= float(nSteps)
	tempAvgAlign /= float(nSteps)
	# Calculating the residual (RMSD) between the avgAlign and tempAvgAlign arrays; tempAvgAlign will become the new avgAlign array; tempAvgCoord will become the new avgCoord 
        residual = RMSD(avgAlign, tempAvgAlign, u_align_atoms)
	rmsd_all = RMSD(avgCoord, tempAvgCoord, u_important_atoms)
	avgCoord = tempAvgCoord
	avgAlign = tempAvgAlign
	ffprint('Steps: %d, RMSD btw alignment atoms: %e, RMSD btw all atoms: %e' %(iteration,residual,rmsd_all))
        iteration += 1

ffprint('Average structure has converged')				# Now have the iteratively aligned avgCoord array, as well as the iteratively aligned (COG-corrected and rotated) allCoord array

# ----------------------------------------
# INITIATE THE PDB FILE TO BE USED TO SAVE THE AVERAGE STRUCTURE
avg_pdb = MDAnalysis.Universe(parameters['pdb'])
avg_important = avg_pdb.select_atoms(parameters['important'])
avg_align = avg_pdb.select_atoms(parameters['alignment'])

# ----------------------------------------
# Print out pdb of average structure
ffprint('Writing a pdb of the average structure.')
avg_important.positions = avgCoord
R,rmsd = rotation_matrix(avg_align.positions,u_align_pos)
avg_important.rotate(R)
avg_important.write('%03d.%03d.avg_structure.pdb' %(parameters['start'],parameters['end']))

ffprint('Finished writing pdb of the average structure')

# ----------------------------------------
# APPENDING INFORMATION TO A SUMMARY FILE
out = open(parameters['out_file'],'a')
out.write('%d   %d   %d\n' %(parameters['start'],parameters['end'], nSteps))
out.close()

