#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

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

input_file = sys.argv[1]
pdb_file = sys.argv[2]

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'

zeros = np.zeros
dot_prod = np.dot
sqrt = np.sqrt
flush = sys.stdout.flush

thresh = 1E-5
maxIter = 100

# ----------------------------------------
# SUBROUTINES:

def ffprint(string):
	print '%s' %(string)
	flush()

# ----------------------------------------
# MAIN PROGRAM:

# LOAD IN AVERAGE SUMMARY FILE AND CALCULATE THE IMPORTANT NUMBERS
average_list = np.loadtxt('%s' %(input_file))

nAverages = len(average_list)
nSteps = np.sum(average_list[:,2])

# LOAD IN THE PDB FILE TO BE USED AS THE AVERAGE STRUCTURE
u = MDAnalysis.Universe(pdb_file)
u_all = u.select_atoms('all')
u_align = u.select_atoms(alignment)
u_important = u.select_atoms('protein or nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')
u_substrate = u.select_atoms('nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')
pos0 = u_align.positions

# GRABBING IMPORTANT NUMBERS FROM THE UNIVERSE
u_important_atoms = len(u_important.atoms)
u_align_atoms = len(u_align.atoms)
u_substrate_res = len(u_substrate.residues)

# ARRAY DECLARATION
all_coord = zeros((nAverages,u_important_atoms,3),dtype=np.float32)
avgCoord = zeros((u_important_atoms,3),dtype=np.float32)
all_align = zeros((nAverages,u_align_atoms,3),dtype=np.float32)
avgAlign = zeros((u_align_atoms,3),dtype=np.float32)

# AVERAGE STRUCTURE ANALYSIS
ffprint('Beginning the averaging of averages process')
for i in range(nAverages):
	ffprint('Loading in the average structure of Trajectories %d to %d' %(average_list[i][0],average_list[i][1]))
	temp = MDAnalysis.Universe('%03d.%03d.avg_structure.pdb' %(average_list[i][0],average_list[i][1]))
	# INITIATING ATOM SELECTIONS
	temp_all = temp.select_atoms('all')
	temp_align = temp.select_atoms(alignment)
	temp_important = temp.select_atoms('protein or nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')
	temp_substrate = temp.select_atoms('nucleic or resname A5 or resname A3 or resname U5 or resname atp or resname adp or resname PHX or resname MG')

	# TRANSLATING AVERAGE STRUCTURES TO THE PDB STRUCTURE
	temp_all.translate(-temp_align.center_of_mass())

	# UNWEIGHTING THE AVERAGE STRUCTURE
	all_coord[i] = temp_important.positions
	all_align[i] = temp_align.positions
	avgCoord += temp_important.positions*average_list[i][2]
	avgAlign += temp_align.positions*average_list[i][2]

avgCoord /= float(nSteps)
avgAlign /= float(nSteps)
ffprint('Finished collecting the average coordinates and averaging the averages...')

# Calculating and Aligning to the average positions
iteration = 0
residual = thresh + 10.0 					# arbitrary assignment greater than thresh
ffprint('Beginning iterative process of calculating average positions and aligning to the average')
while residual > thresh and iteration < maxIter:		
	tempAvgCoord = zeros((u_important_atoms,3),dtype=np.float32)		# zeroing out the tempAvgCoord array every iteration
	tempAvgAlign = zeros((u_align_atoms,3),dtype=np.float32)
	for i in range(nAverages):
		R, d = rotation_matrix(all_align[i,:,:],avgAlign)
		all_align[i,:,:] = dot_prod(all_align[i,:,:],R.T)
		all_coord[i,:,:] = dot_prod(all_coord[i,:,:],R.T)
		tempAvgAlign += all_align[i,:,:]
		tempAvgCoord += all_coord[i,:,:]
	tempAvgCoord /= float(nAverages)
	tempAvgAlign /= float(nAverages)
	residual = RMSD(avgAlign,tempAvgAlign,u_align_atoms)
	rmsd_all = RMSD(avgCoord,tempAvgCoord,u_important_atoms)
	iteration += 1
	avgCoord = tempAvgCoord
	avgAlign = tempAvgAlign
	ffprint('Steps: %d, RMSD btw alignment atoms: %e, RMSD btw all atoms: %e' %(iteration,residual,rmsd_all))

ffprint('Average structure has converged')

u_important.positions = avgCoord
u_important.write('%03d.%03d.average_structure.pdb' %(average_list[0][0],average_list[-1][1]))
ffprint('Finished writing pdb of the average structure')

