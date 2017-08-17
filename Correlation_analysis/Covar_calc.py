# ----------------------------------------
# USAGE:
# ----------------------------------------

# python res_res_correlation.py config_file

# ----------------------------------------
# PREAMBLE:
# ----------------------------------------

import sys
import numpy as np
from numpy.linalg import *
import MDAnalysis
from MDAnalysis.analysis.align import *

zeros = np.zeros
sqrt = np.sqrt
sum = np.sum
eigen = np.linalg.eig
flush = sys.stdout.flush

# ----------------------------------------
# VARIABLE DECLARATION
# ----------------------------------------

config_file = sys.argv[1]

necessary_parameters = ['pdb_file','traj_loc','start','end','average_pdb']
all_parameters = ['pdb_file','traj_loc','start','end','average_pdb','alignment','covar_selection','coarseness','fine_grain_selection','cartesian_correlation_filename','cartesian_average_filename','cartesian_variance_filename','distance_correlation_filename','distance_variance_filename','functionalize_distance_correlation_bool','functionalized_distance_correlation_filename','PCA_bool','PCA_eigenvalues_filename','PCA_eigenvectors_filename','summary_bool','summary_filename']

# ----------------------------------------
# SUBROUTINES:
# ----------------------------------------

def ffprint(string):
	print '%s' %(string)
	flush()

def config_parser(config_file):	# Function to take config file and create/fill the parameter dictionary 
	for i in range(len(necessary_parameters)):
		parameters[necessary_parameters[i]] = ''
	
	# SETTING DEFAULT PARAMETERS FOR OPTIONAL PARAMETERS:
	parameters['alignment'] = 'protein'
	parameters['covar_selection'] = 'protein'
	parameters['coarseness'] = 'COM'
	parameters['fine_grain_selection'] = None
	parameters['cartesian_correlation_filename'] = 'cartesian_correlation.dat'
	parameters['cartesian_average_filename'] = 'cartesian_average.dat'
	parameters['cartesian_variance_filename'] = 'cartesian_variance.dat'
	parameters['distance_correlation_filename'] = 'distance_correlation.dat'
	parameters['distance_variance_filename'] = 'distance_variance.dat'
	parameters['functionalize_distance_correlation_bool'] = False
	parameters['functionalized_distance_correlation_filename'] = 'functionalized_dist_covar.dat'
	parameters['PCA_bool'] = False
	parameters['PCA_eigenvalues_filename'] = 'PCA_eigenvalues_cartesian_covariance.dat' 
	parameters['PCA_eigenvectors_filename'] = 'PCA_eigenvectors_cartesian_covariance.dat' 
	parameters['summary_bool'] = True
	parameters['summary_filename'] = 'water_retention_analysis.summary'

	# GRABBING PARAMETER VALUES FROM THE CONFIG FILE:
	execfile(config_file,parameters)
	for key, value in parameters.iteritems():
		if value == '':
			print '%s has not been assigned a value. This variable is necessary for the script to run. Please declare this variable within the config file.' %(key)
			sys.exit()

	if parameters['coarseness'] not in ['COM','Atomic']:
		print "coarseness parameter does not match an acceptable value. Viable values are 'COM' and 'Atomic'. Killing job."
		sys.exit()

def summary(filename):
	with open(filename,'w') as W:
		W.write('Using MDAnalysis version: %s\n' %(MDAnalysis.version.__version__))
		W.write('To recreate this analysis, run this line:\n')
		for i in range(len(sys.argv)):
			W.write('%s ' %(sys.argv[i]))
		W.write('\n\nParameters used:\n')
		for i in all_parameters:
			W.write('%s = %s \n' %(i,parameters[i]))
		W.write('\n\n')

# ----------------------------------------
# MAIN:
# ----------------------------------------
# CREATING PARAMETER DICTIONARY

parameters = {}
config_parser(config_file)

start = int(parameters['start'])
end = int(parameters['end'])
# ----------------------------------------
# INITIATE THE AVG STRUCTURE; GRAB THE NECESSARY INFORMATION

ffprint('Initiating the average structure universe')
avg = MDAnalysis.Universe(parameters['average_pdb'])
avg_all = avg.select_atoms('all')
avg_align = avg.select_atoms(parameters['alignment'])
avg_all.translate(-avg_align.center_of_mass())
pos0 = avg_align.positions

# ----------------------------------------
# INITIALIZE THE ANALYSIS UNIVERSE; CREATE THE NECESSARY ATOM SELECTIONS

u = MDAnalysis.Universe(parameters['pdb_file'])
u_all = u.select_atoms('all')
u_align = u.select_atoms(parameters['alignment'])
u_covar = u.select_atoms(parameters['covar_selection'])

# ----------------------------------------
# COARSENESS -- Center Of Mass of Residues in the covar_selection

if parameters['coarseness'] == 'COM':
	ffprint('Performing a correlation analysis of the cartesian coordinates of the center of mass of residues defined in the covar_selection parameter.')

	nNodes = u_covar.n_residues
	nDims = nNodes*3			# each node is a point in cartesian space; therefore, the number of dimensions to collect data of is nNodes*3
	if nNodes != avg.select_atoms(parameters['covar_selection']).n_residues:
		ffprint('The number of residues to be used in the covar_selection do not match between the average and analysis universes. Killing job.')
		sys.exit()

	# ----------------------------------------
	# MEMORY DECLARATION
	covariance_array = zeros((nDims,nDims),dtype=np.float64)
	average_array = zeros((nDims),dtype=np.float64)
	variance_array = zeros((nDims),dtype=np.float64)

	# ----------------------------------------
	# TRAJECTORY ANALYSIS
	ffprint('Beginning trajectory analysis.')
	nSteps = 0
	while start <= end:
		ffprint('Loading trajectory %s' %(start))
		u.load_new(parameters['traj_loc'] %(start,start))
		nSteps += len(u.trajectory)
		for ts in u.trajectory:
			u_all.translate(-u_align.center_of_mass())
			R,d = rotation_matrix(u_align.positions,pos0)		# MDAnalysis.analysis.align function
			u_all.rotate(R)
			
			com_array = zeros((nDims),dtype=np.float64)
			for i in range(nNodes):
				j = i*3
				com_array[j:j+3] = u_covar.residues[i].center_of_mass() 
			
			for i in range(nDims):	
				average_array[i] += com_array[i]				# summing over nSteps; x_{i}(t)
				variance_array[i] += com_array[i]**2				# summing over nSteps; x_{i}(t)**2
				for j in range(i,nDims):
					covariance_array[i,j] += com_array[i]*com_array[j]	# summing over nSteps; x_{i}(t) * x_{j}(t)
		start += 1

# ----------------------------------------
# COARSENESS -- Atomic positions of fine_grain_selection in the covar_selection

elif parameters['coarseness'] == 'Atomic':
	ffprint('Performing a covariance analysis of the cartesian coordinates of the fine_grain_selection for the covar_selection.')

	u_fine_grain = u_covar.select_atoms(parameters['fine_grain_selection'])
	nNodes = u_fine_grain.n_atoms
	nDims = nNodes*3			# each node is a point in cartesian space; therefore, the number of dimensions to collect data of is nNodes*3

	avg_covar = avg.select_atoms(parameters['covar_selection'])
	if nNodes != avg_covar.select_atoms(parameters['fine_grain_selection']).n_atoms:
		ffprint('The number of atoms to be analyzed in the fine_grain_selection of the covar_selection do not match between the average and analysis universes. Killing job.')
		sys.exit()

	# ----------------------------------------
	# MEMORY DECLARATION
	covariance_array = zeros((nDims,nDims),dtype=np.float64)
	average_array = zeros((nDims),dtype=np.float64)
	variance_array = zeros((nDims),dtype=np.float64)

	# ----------------------------------------
	# TRAJECTORY ANALYSIS
	ffprint('Beginning trajectory analysis.')
	nSteps = 0
	start = int(parameters['start'])
	while start <= end:
		ffprint('Loading trajectory %s' %(start))
		u.load_new(parameters['traj_loc'] %(start,start))
		nSteps += len(u.trajectory)
		for ts in u.trajectory:
			u_all.translate(-u_align.center_of_mass())
			R,d = rotation_matrix(u_align.positions,pos0)		# MDAnalysis.analysis.align function
			u_all.rotate(R)
			
			pos_array = zeros((nDims),dtype=np.float64)
			for i in range(nNodes):
				j = i*3
				pos_array[j:j+3] = u_covar.atoms[i].position 
			
			for i in range(nDims):
				average_array[i] += pos_array[i]				# summing over nSteps; x_{i}(t)
				variance_array[i] += pos_array[i]**2				# summing over nSteps; x_{i}(t)**2
				for j in range(i,nDims):
					covariance_array[i,j] += pos_array[i]*pos_array[j]	# summing over nSteps; x_{i}(t) * x_{j}(t)
		start += 1


# ----------------------------------------
# FINISHING AVERAGES OF THE AVERAGE, VARIANCE, AND COVARIANCE ARRAYS

covariance_array /= nSteps			# finishing the average over nSteps; <x_{i}(t) * x_{j}(t)>
average_array /= nSteps				# finishing the average over nSteps; <x_{i}(t)>
variance_array /= nSteps			# finishing the average over nSteps; <x_{i}(t)**2>
variance_array -= average_array**2		# finishing the variance analysis; <x_{i}(t)**2> - <x_{i}(t)>**2

# ----------------------------------------
# FINISHING THE CARTESIAN CORRELATION MATRIX OF RESIDUE-RESIDUE PAIRS 

cart_correlation_matrix = zeros((nDims,nDims),dtype=np.float64)
for i in range(nDims):
	for j in range(i,nDims):		# loops through all necessary elements
		covariance_array[i,j] -= average_array[i]*average_array[j]	# finishing the covariance analysis; <x_{i}(t) * x_{j}(t)> - <x_{i}(t)>*<x_{j}(t)>; storing this array for the PCA analysis later on
		cart_correlation_matrix[i,j] = covariance_array[i,j]/sqrt(variance_array[i]*variance_array[j])	    # finishing the correlation analysis; <x_{i}(t) * x_{j}(t)> - <x_{i}(t)>*<x_{j}(t)>/ sqrt((<x_{i}(t)**2> - <x_{i}(t)>**2)(<x_{j}(t)**2> - <x_{j}(t)>**2))
		cart_correlation_matrix[j,i] = cart_correlation_matrix[i,j]					# filling in the bottom triangle of this matrix

# ----------------------------------------
# OUTPUTING CARTESIAN DATA
with open(parameters['cartesian_correlation_filename'],'w') as f:
	np.savetxt(f,cart_correlation_matrix)

with open(parameters['cartesian_average_filename'],'w') as f:
	np.savetxt(f,average_array)

with open(parameters['cartesian_variance_filename'],'w') as f:
	np.savetxt(f,variance_array)

# ----------------------------------------
# FINISHING THE DISTANCE CORRELATION MATRIX OF RESIDUE-RESIDUE PAIRS
distance_correlation_matrix = zeros((nNodes,nNodes),dtype=np.float64)
distance_variance_array = zeros((nNodes),dtype=np.float64)
for i in range(nNodes):
	dim1 = i*3
	distance_variance_array[i] = sum(variance_array[dim1:dim1+3])						# summing the variances of the dimensions of node i; <r_{i}(t)**2>
	for j in range(i,nNodes):
		dim2 = j*3
		distance_correlation_matrix[i,j] = covariance_array[dim1,dim2] + covariance_array[dim1+1,dim2+1] + covariance_array[dim1+2,dim2+2]	# taking the trace of the covariance analysis matrix for the desired dimensions of nodes i and j; <r_{i}(t) dot r_{j}(t)> - <r_{i}(t)> dot <r_{j}(t)> 

for i in range(nNodes):
	for j in range(i,nNodes):
		distance_correlation_matrix[i,j] /= sqrt(distance_variance_array[i]*distance_variance_array[j])	    # finishing the coorelation analysis; <r_{i}(t) dot r_{j}(t)> - <r_{i}(t)> dot <r_{j}(t)>/  sqrt((<r_{i}(t)**2> - <r_{i}(t)>**2)(<r_{j}(t)**2> - <r_{j}(t)>**2))
		distance_correlation_matrix[j,i] = distance_correlation_matrix[i,j]	# filling in the bottom triangle of this matrix

# ----------------------------------------
# OUTPUTING DISTANCE DATA
with open(parameters['distance_correlation_filename'],'w') as f:
	np.savetxt(f,distance_correlation_matrix)

with open(parameters['distance_variance_filename'],'w') as f:
	np.savetxt(f,distance_variance_array)

# ----------------------------------------
# FUNCTIONALIZING THE DISTANCE CORRELATION MATRIX (FOR USE IN WISP AND VISUALIZATION)
if parameters['functionalize_distance_correlation_bool']:
	ffprint('Beginning to functionalize the distance covar matrix.')
	for node1 in range(nNodes):
		for node2 in range(node1,nNodes):
			distance_correlation_matrix[node1,node2] = -np.log(np.fabs(distance_correlation_matrix[node1,node2]))
			distance_correlation_matrix[node2,node1] = distance_correlation_matrix[node1,node2]
	
	# OUTPUTING THE DIST COVAR ARRAY
	with open(parameters['functionalized_distance_correlation_filename'],'w') as f:
		np.savetxt(f,distance_correlation_matrix)
else: 
	ffprint('Functionalize != True. Not functionalizing (taking -log(|C_ij|)) the distance correlation matrix.')

# ----------------------------------------
# PCA ANALYSIS OF CARTESIAN COVAR ARRAY -- Be Warned! This is time consuming
if parameters['PCA_bool']:
	ffprint('Beginning PCA analysis of the Cartesian covariance matrix.')
	eigval,eigvec = eigen(covariance_array)
	idx = eigval.argsort()[::-1]
	eigval = eigval[idx]
	
	nVec = len(eigvec)
	cumulative_eigval = zeros(nVec,dtype=np.float64)
	total_eigval = 0
	for i in range(nVec):
		total_eigval += eigval[i]
		cumulative_eigval[i] = total_eigval
	
	with open(parameters['PCA_eigenvalues_filename'],'w') as f:
		for i in range(nVec):
			f.write('%f   %f   %f   %f\n' %(eigval[i],eigval[i]/total_eigval,cumulative_eigval[i],cumulative_eigval[i]/total_eigval))
	
	with open(parameters['PCA_eigenvectors_filename'],'w') as f:
		for i in range(nVec):
			for j in range(nVec):
				f.write('%f   ' %(eigvec[j,i]))		# Writing each vector on one row/line now, instead of the vectors corresponding to columns in the eigvec array...; NOT projecting covar array onto the eigenvectors (do so outside of this damn script)
			f.write('\n')
else:
	ffprint('PCA_bool != True. Not performing a PCA analysis on the Cartesian Covar matrix.')

if parameters['summary_bool']:
	summary(parameters['summary_filename'])

