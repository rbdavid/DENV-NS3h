#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# ----------------------------------------
# USAGE:

# ./rmsd.ref.py pdb_file trajectory_location start end ref_pdb_file system_descriptor

# ----------------------------------------
# PREAMBLE:

import numpy as np
import sys
import os
import MDAnalysis
from MDAnalysis.analysis.align import *
from sel_list import *
from distance_functions import *

pdb = sys.argv[1]			# point to a pdb or prmtop or psf file (untested for both prmtop and psf files)
traj_loc = sys.argv[2]			# point to the location of the trajectory files
start = int(sys.argv[3])		# integer describing 
end = int(sys.argv[4])
ref_pdb = sys.argv[5]
system = sys.argv[6]

alignment = 'protein and name CA and (resid 20:25 or resid 50:55 or resid 73:75 or resid 90:94 or resid 112:116 or resid 142:147 or resid 165:169 or resid 190:194 or resid 214:218 or resid 236:240 or resid 253:258 or resid 303:307)'

nSel = len(sel)

flush = sys.stdout.flush

# ----------------------------------------
# FUNCTIONS:

def ffprint(string):
	print '%s' %(string)
        flush()

def summary(nSteps):
	sum_file = open('%s.rmsd.summary' %(system),'w')
	sum_file.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
	sum_file.write('To recreate this analysis, run this line:\n')
	for i in range(len(sys.argv)):
		sum_file.write('%s ' %(sys.argv[i]))
	sum_file.write('\n\n')
	sum_file.write('output is written to:\n')
	sum_file.write('	%s.rmsd.dat\n' %(system))
	sum_file.write('\nNumber of steps analyzed: %d\n' %(nSteps))
	for i in range(nSel):
		sum_file.write('	%02d   %s   %s\n' %(i,sel[i][0],sel[i][1]))
	sum_file.close()

# ----------------------------------------
# MAIN:

ref = MDAnalysis.Universe('%s' %(ref_pdb))
ref_align = ref.select_atoms(alignment)
ref_all = ref.select_atoms('all')
ref_all.translate(-ref_align.center_of_mass())
pos0 = ref_align.positions

u = MDAnalysis.Universe('%s' %(pdb))
u_align = u.select_atoms(alignment)
u_all = u.select_atoms('all')

if len(u_align) != len(ref_align):
	ffprint('Alignment atom selections do not have the same number of atoms.')
	sys.exit()

rest = u.select_atoms('not (resname WAT or resname Na+ or resname Cl- or protein)')
num_res = len(rest.residues)

# make selections to compute RMSD for
u_sel = []
ref_sel = ['']*nSel
for i in range(nSel):
	selection = sel[i][1]

	ref_temp = ref.select_atoms(selection)
	u_temp = u.select_atoms(selection)

	u_sel.append([int(u_temp.n_atoms),u_temp])
	
	if len(u_sel[i][1]) != len(ref_temp):
		ffprint('Number of atoms do not match for selection %02d' %(i))
		sys.exit()

	ref_sel[i] = ref_temp.positions

# open output files
out1 = open('%s.rmsd.dat' %(system), 'w')
nSteps = 0

ffprint('Beginning trajectory analysis')
while start <= end:
	ffprint('Loading trajectory %s' %(start))
	u.load_new('%sproduction.%s/production.%s.dcd' %(traj_loc,start,start))
	nSteps += len(u.trajectory)
	
	# Loop through trajectory
	for ts in u.trajectory:
		# obtain dimension values to be used for unwrapping atoms
		dims = u.dimensions[:3]
	
		# Align to reference (moves COM of alignment to origin)
		u_all.translate(-u_align.center_of_mass())

		# Fix the wrapping issues
		for i in range(num_res):
			COM = np.zeros(3)
			res = rest.residues[i]
			# Calculate the COM of residues;
			COM = res.center_of_mass()
			# CALCULATING AND APPLYING THE TRANSLATIONAL MATRIX TO RESIDUE i
			t = wrapping(COM,dims)
			res.atoms.translate(t)

		# Calculate the rotational matrix to align u to the ref
		R, rmsd = rotation_matrix(u_align.positions, pos0)
		# Apply rotation matrix to atoms within u
		u_all.rotate(R)

		# loop through selections and compute RMSD
		for i in range(nSel):
			u_coords = u_sel[i][1].positions
			ref_coords = ref_sel[i]
			
			rmsd = RMSD(u_coords,ref_coords,u_sel[i][0])
			out1.write('%10.6f   ' %(rmsd))

		out1.write('\n')
	start += 1

out1.close()
summary(nSteps)

