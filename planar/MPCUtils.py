import numpy as np
import scipy.sparse as sparse # for testing
import osqp
import sys

'''
NOTE: Dimensions ---
Ad = nx*nx, Bd = nx*nu
Ax = (N+1)*nx * (N+1)*nx
Bu = (N+1)*nx * N*nu
Aeq = (N+1)*nx * ((N+1)*nx+N*nu)
x = [y0, y1, ..., yN, u1, ..., uN] of size (N+1)*nx + N*nu
After solving, apply u1 (MPC)
# print(P.shape)
# print(q.shape)
# print(A.shape)
# print(l.shape)
# print(leq.shape)
# print(lineq.shape)
'''

class MPCHelper:
	def __init__(self, model, N, wx, wu, **settings):
		# Pass in nx=#states, nu=#inputs, and N=prediction horizon
		# model must have:
		# nx, nu, N
		# getLimits()
		# getLin()

		self.m = model
		self.N = N

		# Constraints
		umin, umax, xmin, xmax = model.getLimits()

		# Create an OSQP object
		self.prob = osqp.OSQP()
		
		# Objective function
		self.Q = sparse.diags(wx)
		self.QN = self.Q
		R = sparse.diags(wu)

		# Cast MPC problem to a QP: x = (x(0),x(1),...,x(N),u(0),...,u(N-1))
		# - quadratic objective
		P = sparse.block_diag([sparse.kron(sparse.eye(self.N), self.Q), self.QN, sparse.kron(sparse.eye(N), R)]).tocsc()
		# - input and state constraints
		Aineq = sparse.eye((N+1)*self.m.nx + N*self.m.nu)
		lineq = np.hstack([np.kron(np.ones(N+1), xmin), np.kron(np.ones(N), umin)])
		uineq = np.hstack([np.kron(np.ones(N+1), xmax), np.kron(np.ones(N), umax)])

		# Create the CSC A matrix manually. 
		conA = CondensedA(self.m.nx, self.m.nu, N)
		self.A = sparse.csc_matrix((conA.data, conA.indices, conA.indptr))
		# make up xr for now - will get updated
		xr = np.zeros(self.m.nx)
		x0 = np.zeros_like(xr)
		q = np.hstack([np.kron(np.ones(N), -self.Q.dot(xr)), -self.QN.dot(xr), np.zeros(N*self.m.nu)])
		# Initial state will get updated
		leq = np.hstack([-x0, np.zeros(self.N*self.m.nx)])
		ueq = leq
		self.l = np.hstack([leq, lineq])
		self.u = np.hstack([ueq, uineq])
		
		# Full so that it is not made sparse. prob.update() cannot change the sparsity structure
		# Setup workspace
		self.prob.setup(P, q, self.A, self.l, self.u, warm_start=True, **settings)#, eps_abs=1e-05, eps_rel=1e-05
		self.ctrl = np.zeros(self.m.nu)

	def update(self, x0, xr, dt):
		# Pass xr = goal

		# - linear objective
		q = np.hstack([np.kron(np.ones(self.N), -self.Q.dot(xr)), -self.QN.dot(xr), np.zeros(self.N*self.m.nu)])
		
		# -- Updating A --
		
		for ti in range(self.N):
			if len(x0.shape) > 1:
				# a whole trajectory has been provided
				Ad, Bd = self.m.getLin(x0[ti,:], self.ctrl, dt)
			elif ti == 0:
				# only the current state provided; only need to call once
				Ad, Bd = self.m.getLin(x0, self.ctrl, dt)
			# Update the LTV dynamics
			cscUpdateDynamics(self.A, self.N, ti, Ad=Ad, Bd=Bd)
		# Update initial state
		if len(x0.shape) > 1:
			self.l[:self.m.nx] = -x0[0,:]
			self.u[:self.m.nx] = -x0[0,:]
		else:
			self.l[:self.m.nx] = -x0
			self.u[:self.m.nx] = -x0
		self.prob.update(l=self.l, u=self.u, q=q, Ax=self.A.data)
		# print(A.data.shape)

		# Solve
		res = self.prob.solve()

		# Check solver status
		if res.info.status != 'solved':
			raise ValueError('OSQP did not solve the problem!')

		# Apply first control input to the plant
		self.ctrl = res.x[-self.N*self.m.nu:-(self.N-1)*self.m.nu]

		# if self.ctrl[0] < -1e-6:

		# 	# FIXME: u seems out of bounds
		# 	Ax = self.A @ res.x
		# 	print((Ax - self.l)[-self.N*self.m.nu:-(self.N-1)*self.m.nu])
		# 	print((self.u - Ax)[-self.N*self.m.nu:-(self.N-1)*self.m.nu])
		# 	print("HELLO")
		# 	sys.exit(0)

		return self.ctrl


'''
These are all for the condensed A matrix
CSC matrix
three NumPy arrays: indices, indptr, data
In OSQP: sparse structure cannot be changed, i.e. can only change data
Ax_new: Vector of new elements in A->x
Ax_new_idx: Index mapping new elements to positions in A->x
A->x must be the data vector
Assuming that indices, indptr can't be changed
'''

class CondensedA:

	def __init__(self, nx, nu, N, val = 0):
		# Pass in nx=#states, nu=#inputs, and N=prediction horizon
		self.nx = nx
		self.nu = nu
		self.N = N

		# create the indptr and index arrays for sparse CSC format
		self.data = np.zeros((self.nnz))
		self.indices = np.zeros_like(self.data, dtype=np.int32)
		numcols = self.shape[1]
		self.indptr = np.zeros((numcols + 1), dtype=np.int32)

		# NOTE: in the code below, i refers to row index, j to col index
		# similar for bj for block col index
		
		# populate indptr
		si = 0
		for j in range(numcols):
			if j < self.N * self.nx:
				si += nx + 2
			elif j < (self.N + 1) * self.nx:
				si += 2
			else:
				si += self.nx + 1
			self.indptr[j+1] = si
		
		# populate row indices
		self._offs = -1 # internal val to count through the data vector
		for j in range(self.shape[1]):
			# In each column, the Aeq blocks come first
			if j < (self.N + 1) * self.nx:
				# Leftmost cols
				bj = int(np.floor(j / self.nx)) # block col index - goes from [0,N+1)
				# the block diagonal -I
				self.addNZ(j, j, -1)
					# The sub-block-diagonal Ad under the -I -- only in the left column block
				if j < self.N * self.nx:
					for i in range(self.nx * (bj + 1), self.nx * (bj + 2)):
						self.addNZ(j, i, val)
			else:
				# Right block column
				bj = int(np.floor((j - (self.N + 1) * self.nx) / self.nu))# block col index for the right - goes from [0,N)
				# The Bd blocks
				for i in range(self.nx * (bj + 1), self.nx * (bj + 2)):
					self.addNZ(j, i, val)

			# The identity in Aineq
			self.addNZ(j, (self.N + 1) * self.nx + j, 1)


			# offs = self.matrixIdxToOffset(i, j)
			# if offs >= 0:
			# 	self.indices[offs] = j

		# sparse structure is now fixed, but data can be updated
	
	def addNZ(self, j, i, val):
		self._offs += 1
		self.indices[self._offs] = i # row index
		self.data[self._offs] = val

	# def blockToMatrixIdx(self, ib, jb, i, j):
	# 	# Return matrix indices from the block index (ib, jb) and element index within the block (i, j)
	# 	if ib < 2*(N+1):
	# 		rowi = self.nx * ib + i
	# 	else:
	# 		rowi = 2 * (self.N+1) * self.nx + (ib - 2 * (self.N+1)) * self.nu + i

	# 	if jb < (N+1):
	# 		rowi = self.nx * jb + j
	# 	else:
	# 		rowi = (self.N+1) * self.nx + (jb - (self.N+1)) * self.nu + j
		
	# 	return rowi, coli

	def get_nnz(self):
		# For testing
		return self.N * self.nx * (self.nx + 2) + self.nx * 2 + self.N * self.nu * (self.nx + 1)
	def get_shape(self):
		# For testing
		return (2 * (self.N + 1) * self.nx + self.N * self.nu), (self.N + 1) * self.nx + self.N * self.nu

	nnz = property(get_nnz, None, None, "Number of nonzero entries")
	shape = property(get_shape, None, None, "Number of nonzero entries")

# These work for sparse as well as CondensedA

def cscUpdateElem(obj, i, j, val):
	# will update an element; do nothing if that entry is zero
	# works on scipy sparse as well as CondensedA

	# indptr has #cols elements, and the entry is the index into data/indices for the first element of that column
	offs = obj.indptr[j]
	# 
	if j < obj.shape[1] - 1:
		maxOffsThisCol = obj.indptr[j+1]
	else:
		maxOffsThisCol = len(obj.indices)
	while offs < maxOffsThisCol:
		if obj.indices[offs] == i:
			obj.data[offs] = val
			return
		offs += 1

def cscUpdateDynamics(obj, N, ti, Ad=None, Bd=None):
	# pass a block to update, and a matrix to go in that block, and it will return a tuple (data, indices) of the same length that can go in to update a sparse matrix
	
	if Ad is not None:
		nx = Ad.shape[0]
		for i in range(nx):
			for j in range(nx):
				cscUpdateElem(obj, nx * (ti + 1) + i, nx * ti + j, Ad[i,j])

	if Bd is not None:
		nx, nu = Bd.shape
		for i in range(nx):
			for j in range(nu):
				cscUpdateElem(obj, nx * (ti + 1) + i, (N + 1) * nx + nu * ti + j, Bd[i,j])
		

if __name__ == "__main__":
	print("Testing CondensedA...")

	#
	nx = 3
	nu = 2
	N = 4
	# create instance
	conA = CondensedA(nx, nu, N, val=1)
	# wider display
	# np.set_printoptions(edgeitems=30, linewidth=100000)
	
	# Create dense and then use scipy.sparse
	Ad = np.ones((nx, nx))
	Bd = np.ones((nx, nu))
	# Ad, Bd = getLin(x0, ctrl, dt)
	Ax = sparse.kron(sparse.eye(N+1),-sparse.eye(nx)) + sparse.kron(sparse.eye(N+1, k=-1), Ad)
	Bu = sparse.kron(sparse.vstack([sparse.csc_matrix((1, N)), sparse.eye(N)]), Bd)
	Aeq = sparse.hstack([Ax, Bu])
	Aineq = sparse.eye((N+1)*nx + N*nu)
	A = sparse.vstack([Aeq, Aineq]).tocsc()

	# Tests ---

	assert(conA.shape == A.shape)
	assert(conA.nnz == A.nnz)
	assert((conA.indptr == A.indptr).all())
	assert((conA.indices == A.indices).all())
	# both were filled with ones
	assert((conA.data == A.data).all())

	#  Usage: update blocks on a sparse, and a CondensedA, and compare to a newly created sparse
	Ad2 = np.full((nx, nx), 123)
	Bd2 = np.full((nx, nu), -456)
	Ax = sparse.kron(sparse.eye(N+1),-sparse.eye(nx)) + sparse.kron(sparse.eye(N+1, k=-1), Ad2)
	Bu = sparse.kron(sparse.vstack([sparse.csc_matrix((1, N)), sparse.eye(N)]), Bd2)
	Aeq = sparse.hstack([Ax, Bu])
	A2 = sparse.vstack([Aeq, Aineq]).tocsc()
	for ti in range(N):
		cscUpdateDynamics(A, N, ti, Ad=Ad2, Bd=Bd2)
		cscUpdateDynamics(conA, N, ti, Ad=Ad2, Bd=Bd2)
	# test update
	assert((conA.data == A2.data).all())
	assert((A.data == A2.data).all())

	# print(conA.indices)
	# print(A)
	# print(conA.data)
	# cscUpdateElem(conA, 6, 0, 123)
	# print(conA.data)

	# print(A.data)
	# print(conA.shape)
	# print(conA.indptr.shape)
	# # cscUpdateDynamics(A, N, 9, Bd=np.full((nx, nu), 123))
	# print(A.toarray())

	print("All passed")
