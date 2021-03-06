# ----------------------------------------
# USAGE:

# python water_positioning.py config_file 

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import *
from distance_functions import *

flush = sys.stdout.flush
zeros = np.zeros
sqrt = np.sqrt
sums = np.sum
square = np.square
dot = np.dot
arccos = np.arccos
mean = np.mean
unit_conversion = 180./np.pi

config_file = sys.argv[1]

# ----------------------------------------
# FUNCTIONS: 

def ffprint(string):
	print '%s' %(string)
        flush()

necessary_parameters = ['pdb','traj_loc','start','end','pocket_selection','pocket_radius','wat_resname','wat_O_name','substrate_atom1','substrate_atom2']
all_parameters = ['pdb','traj_loc','start','end','pocket_selection','pocket_radius','wat_resname','wat_O_name','substrate_atom1','substrate_atom2','Wrapped','write_summary','summary_filename','nucl_wat_outputname','avg_wat_outputname','center_of_geometry_filename'] 

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''

	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['write_summary'] = False
	parameters['summary_filename'] = None
	parameters['nucl_wat_outputname'] = 'nucleophilic_waters.dat'
	parameters['avg_wat_outputname'] = 'average_waters.dat'
	parameters['center_of_geometry'] = 'COG.xyz'

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

# ----------------------------------------
# MAIN:
# CREATING PARAMETER DICTIONARY
parameters = {}
config_parser(config_file)

# ----------------------------------------
# LOAD IN THE ANALYSIS UNIVERSE AND CREATE THE NECESSARY ATOM SELECTIONS
ffprint('Loading Analysis Universe.')
u = MDAnalysis.Universe(parameters['pdb'])
u_all = u.select_atoms('all')
wat = u.select_atoms(parameters['wat_resname'])
u_pocket = u.select_atoms(parameters['pocket_selection'])
atom1 = u.select_atoms(parameters['substrate_atom1'])
atom2 = u.select_atoms(parameters['substrate_atom2'])

# ----------------------------------------
# ANALYSIS OF TRAJECTORIES
start = int(parameters['start'])
end = int(parameters['end'])
timestep = 0
with open(parameters['nucl_wat_outputname'],'w') as W, open(parameters['avg_wat_outputname'],'w') as X, open(parameters['center_of_geometry_filename'],'w') as Y:
	ffprint('Beginning trajectory analysis')
	while start <= end:
		ffprint('Loading trajectory %s' %(start))
		u.load_new(parameters['traj_loc'] %(start,start))
		# Loop through trajectory
		for ts in u.trajectory:
	                # Obtaining COG of pocket; moving origin to this point
			t = u_pocket.center_of_geometry()
			Y.write('1\n  generated by MDAnalysis and RBD\n X         %10.4f         %10.4f         %10.4f\n' %(t[0], t[1], t[2]))	#Writing an xyz trajectory of the center of geometry of the binding pocket; the COG particle is labeled as a dummy atom X
			u_all.translate(-t)

                        # Wrap waters around the center of the NTPase active site if trajectories are not wrapped already
			if not parameters['Wrapped']:
				dims = u.dimensions[:3]	# obtain dimension values to be used for wrapping atoms
				dims2 = dims/2.0
				for i in range(nWats):
					temp = wat.reidues[i].atom[0].position
					t = wrapping(temp,dims,dims2)
					wat.residues[i].translate(t)
		
			pocket_waters = wat.select_atoms('byres point 0 0 0 %d'%(parameters['pocket_radius']))
			nRes = pocket_waters.n_residues
			X.write('%d\n'%(nRes))

			pos1 = atom1.positions[0]
			pos2 = atom2.positions[0]
			dist,dist2 = euclid_dist(pos1,pos2)
			bond_vector = (pos2 - pos1)/dist

			for i in range(nRes):
				res = pocket_waters.residues[i]
				ox_pos = res.select_atoms('name %s'%(parameters['wat_O_name'])).positions[0]
				dist, dist2 = euclid_dist(pos2,ox_pos)
				attack_vector = (pos2 - ox_pos)/dist
				attack_angle = arccos(dot(bond_vector,attack_vector))*unit_conversion
				W.write('%d   %d   %f   %f\n'%(timestep,i,dist,attack_angle))
			
			timestep += 1

		start += 1

if parameters['write_summary']:
	summary()

