import numpy as np
import sys
np.set_printoptions(linewidth=999)


shape_funcs_d = lambda ksi: np.asarray([ -0.5 + ksi, 0.5 + ksi, -2.*ksi])
shape_funcs = lambda ksi: np.asarray([ -0.5*(ksi*(1-ksi)), 0.5*(ksi*(1.+ksi)), 1.-(ksi**2)])

# Element object to store relevant information on node numbers, positions, etc.
class Element:
	def __init__(self, enum, gx1, gx2, gx3, lx1, lx2, n=2):
		assert n == 1 or n == 2
		self.gx1 = gx1
		self.gx2 = gx2
		self.gx3 = gx3
		self.lx1 = lx1
		self.lx2 = lx2
		self.lx3 = (lx1+lx2)/2.
		self.enum = enum
		self.le = abs(lx2-lx1)
		self.weights = [2., 1.]
		self.n = n
		
		if n == 1:
			self.positions = [0.]
		elif n == 2:
			self.positions = [-3**-0.5,3**-0.5]

		self.k = self.get_local_stiffness()	
		self.f = self.get_local_loading()	


	# Argument is to take in number of integration points
	def get_local_stiffness(self):
		k = np.zeros((3,3))
	
		for p in self.positions:
			k += np.outer(shape_funcs_d(p),shape_funcs_d(p))*self.weights[self.n-1]
		
		return 2.*k/self.le
			
	# m represents the size of the global stiffness matrix
	def get_global_stiffness(self, m, DISPLACEMENT_BC=None):
		global_node_list = [self.gx1, self.gx2, self.gx3]
		partial_global = np.zeros((m,m))

		idxs_o = np.dstack(np.meshgrid(global_node_list, global_node_list)).reshape((len(global_node_list)**2,2))
		idxs = idxs_o[~((idxs_o[:,0]==DISPLACEMENT_BC) | (idxs_o[:,1]==DISPLACEMENT_BC))]
		partial_global[idxs[:,0], idxs[:,1]] = self.get_local_stiffness().T.flatten()\
				[~((idxs_o[:,0]==DISPLACEMENT_BC) | (idxs_o[:,1]==DISPLACEMENT_BC))]

		return partial_global


	# m represents the size of the global stiffness matrix
	def get_global_loading(self, m, DISPLACEMENT_BC=None):
		global_node_list = np.array([self.gx1, self.gx2, self.gx3])
		partial_global = np.zeros(m)

		idxs = global_node_list[~(global_node_list==DISPLACEMENT_BC)]
		partial_global[idxs] = self.get_local_loading()[~(global_node_list==DISPLACEMENT_BC)]

		return partial_global



	# Argument is to take in number of integration points
	def get_local_loading(self):
		f = np.zeros(3)
		
		for p in self.positions:
			f += shape_funcs(p)*( (self.le/2.)*p + (self.lx1+self.lx2)/2.)*self.weights[self.n-1]

		return f*self.le/2.


	def get_fe_strains_IP(self, global_d):
		idx = np.array([self.gx1, self.gx2, self.gx3])
		local_d = global_d[idx]
		strains = np.zeros(2)
		for i, p in enumerate(self.positions):
			strains[i] = np.dot(shape_funcs_d(p), local_d) * (2./self.le)
		return strains


	def get_exact_strains_IP(self, global_d):
		strains = np.zeros(2)
		for i, p in enumerate(self.positions):
			x = np.dot(shape_funcs(p), np.array([self.lx1, self.lx2, self.lx3]))
			strains[i] = -0.5*x**2
		return strains

	def get_exact_disp(self):

		disp = np.array([(1./6.)*(-x**3. + 1.) for x in np.array([self.lx1, self.lx2, self.lx3])])

		return disp

##################################################

# We initially start with a 9x9 since we include the node on which displacement is fixed. We will resolve after the matrix
# is fully constructed by removing the corresponding row. It's an easy implementation this way. This is the same logic for 
# the global loading vector as well.
global_stiffness = np.zeros((8,8))
global_stiffness_1 = np.zeros((8,8))
global_loading = np.zeros(8)
DISPLACEMENT_BC = 8

e0 = Element(0,0,1,2,0.,0.25, n=2)
e1 = Element(1,1,3,4,0.25,0.5, n=2)
e2 = Element(2,3,5,6,0.5,0.75, n=2)
e3 = Element(3,5,8,7,0.75,1., n=2)
elem_list = [e0,e1,e2,e3]


e01 = Element(0,0,1,2,0.,0.25, n=1)
e11 = Element(1,1,3,4,0.25,0.5, n=1)
e21 = Element(2,3,5,6,0.5,0.75, n=1)
e31 = Element(3,5,8,7,0.75,1., n=1)
elem_list_1 = [e01,e11,e21,e31]

# Assemble all together
for e in elem_list:
	global_stiffness += e.get_global_stiffness(8, DISPLACEMENT_BC=DISPLACEMENT_BC)
	global_loading += e.get_global_loading(8, DISPLACEMENT_BC=DISPLACEMENT_BC)

for e in elem_list_1:
	global_stiffness_1 += e.get_global_stiffness(8, DISPLACEMENT_BC=DISPLACEMENT_BC)
		

# We now solve the system using Gaussian Elimination

# A fast gaussian elimination solver I came up with it that parallelizes operations for fast runtime
def gaussian(A, b):
	# Fast Forward Elimination
	for k in np.arange(len(A)-1):
		C = A[(k+1):,k] / A[k,k]
		A[(k+1):,k:] -= np.outer(C, A[k,k:])
		b[(k+1):] -= C * b[k]

	# Fast Back Substitution
	b[-1] /= A[-1,-1]
	for k in -1*np.arange(1,len(A)):
		b[k-1] = (b[k-1]-np.dot(A[k-1,k:],b[k:]))/A[k-1,k-1]

	return b

global_displacement = gaussian(global_stiffness, global_loading)

print "\n\nGlobal Stiffness Matrix: \n", global_stiffness
print "\n\nGlobal Stiffness Matrix (n=1): \n", global_stiffness_1
print "\n\nGlobal Loading Vector: \n", global_loading

global_displacement = np.insert(global_displacement,DISPLACEMENT_BC,0.)
exact_displacement = np.zeros(9)
for e in elem_list:
	exact_displacement[[e.gx1, e.gx2, e.gx3]] = e.get_exact_disp()
print '\n\nFEM Displacement: ', global_displacement
print 'Exact Displacement: ', exact_displacement 


print "\n\nFinite Element Strains: \n", 
for i, e in enumerate(elem_list):
	print "\nElement ", i, ":"
	print e.get_fe_strains_IP(global_displacement)

print "\n\nExact Strains: \n",
for i, e in enumerate(elem_list):
	print "\nElement ", i, ":"
	print e.get_exact_strains_IP(global_displacement)


#


