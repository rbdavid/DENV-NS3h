# ----------------------------------------
# USAGE:

# python rmsd.ref.py config_file

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

flush = sys.stdout.flush

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS: 

def ffprint(string):
	print '%s' %(string)
        flush()

necessary_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file']
all_parameters = ['ref_pdb','pdb','traj_loc','start','end','Wrapped','outputname','selection_file','alignment','write_summary','summary_filename','selection_output']
def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment'] = 'protein'
	parameters['write_summary'] = False 
	parameters['summary_filename'] = 'rmsd.summary'
	parameters['summary_filename'] = 'selections.txt'

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
		f.write('output is written to:\n')
		f.write('	%s\n' %(parameters['outputname']))
		f.write('\nNumber of steps analyzed: %d\n' %(nSteps))
		f.write('\nAtom selections analyzed have been written out to %s\n' %(parameters['selection_output']))

# ----------------------------------------
# MAIN:
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOAD IN THE REFERENCE AND ANALYSIS UNIVERSES
ref = MDAnalysis.Universe(parameters['ref_pdb'])
ref_align = ref.select_atoms(parameters['alignment'])
ref_all = ref.select_atoms('all')
ref_all.translate(-ref_align.center_of_mass())
pos0 = ref_align.positions

u = MDAnalysis.Universe(parameters['pdb'])
u_align = u.select_atoms(parameters['alignment'])
u_all = u.select_atoms('all')
if not parameters['Wrapped']:
	rest = u.select_atoms('not (resname WAT Na+ Cl- or protein)')
	rest_nRes = rest.n_residues

if u_align.n_atoms != ref_align.n_atoms:
	ffprint('Alignment atom selections do not have the same number of atoms.')
	sys.exit()

# ----------------------------------------
# READING SELECTION LIST FILE

execfile(parameters['selection_file'])
nSel = len(sel)

# ----------------------------------------
# CREATION OF ATOM SELECTIONS FOR RMSD ANALYSIS

selection_list = []     # saving u universe selections corresponding to selections made in selection_file
nAtoms = []             # saving num. of atoms in u universe selections
ref_pos = []            # saving positions of atoms in ref universe selections
with open(parameters['selection_output'],'w') as f:
	for i in range(nSel):
		u_temp = u.select_atoms(sel[i][1])
		selection_list.append(u_temp)
		nAtoms.append(u_temp.n_atoms)
		ref_temp = ref.select_atoms(sel[i][1])
		ref_pos.append(ref_temp.positions)
		if u_temp.n_atoms != ref_temp.n_atoms:
			ffprint('Number of atoms do not match for selection %02d %s' %(i,sel[i][0]))
			sys.exit()

		f.write("%02d   %s   '%s'\n" %(i,sel[i][0],sel[i][1]))      # outputting details about selections made

count = len(selection_list)
	
nSteps = 0
start = int(parameters['start'])
end = int(parameters['end'])
with open(parameters['outputname'],'w') as f:
	ffprint('Beginning trajectory analysis')
	while start <= end:
		ffprint('Loading trajectory %s' %(start))
		u.load_new(parameters['traj_loc'] %(start,start))
		nSteps += len(u.trajectory)
		# Loop through trajectory
		for ts in u.trajectory:
			# Align to reference (moves COM of alignment to origin)
			u_all.translate(-u_align.center_of_mass())
		
			# CALCULATIONS that are unnecessary if the trajectory is wrapped.
			if not parameters['Wrapped']:		# Test to see if the 'Wrapped' key is equal to False
				dims = u.dimensions[:3]		
				dims2 = dims/2.0
	
				for i in range(rest_nRes):
					COM = rest.residues[i].center_of_mass()
					t = wrapping(COM,dims,dims2)
					rest.residues[i].atoms.translate(t)
	
			# Calculate the rotational matrix to align u to the ref
			R, rmsd = rotation_matrix(u_align.positions, pos0)
			# Apply rotation matrix to all atoms within u
			u_all.rotate(R)
	
			# loop through selections and compute RMSD
			for i in range(count):
				temp_coords = selection_list[i].positions
				
				rmsd = RMSD(temp_coords,ref_pos[i],nAtoms[i])
				f.write('%10.6f   ' %(rmsd))
	
			f.write('\n')
		start += 1

if parameters['write_summary']:
	summary()

