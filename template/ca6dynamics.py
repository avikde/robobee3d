import autograd.numpy as np

# model params
ycp = 10 # mm
mb = 100
ixx = iyy = 3333
izz = 1000
g = 9.81e-3
M = np.diag([mb, mb, mb, ixx, iyy, izz])

# function aeroWrenchAffine(f)
# 	# params?
#     CLmax = 1.8
#     CDmax = 3.4
#     CD0 = 0.4
# 	ρ = 1.225e-3 # [mg/(mm^3)]
# 	Aw = 54.4 #mm^2
# 	cbar = 3.2 #mm
# 	r2h = 0.551

# 	# calculated
# 	AR = Aw/cbar^2
# 	k = 1/2 * ρ * Aw^2 * AR * r2h^2
# 	kt = k * CLmax * π^2
# 	Lw = Aw/cbar
# 	ycp = 0.5 * Lw
	
# 	Ψ = 1.0 # TODO:

# 	Fz = kt * f^2 * cos(Ψ)*sin(Ψ)
# 	return [Fz, ycp*Fz]
# end

def wrenchMap(u):
    u1L, u2L, u3L, u1R, u2R, u3R = u # unpack
    return np.array([u3L + u3R, 0.0, u1L + u1R, 
        (u1L - u1R) * ycp, -u1L*u2L - u1R*u2R, (-u3L + u3R)*ycp])

def dynamicsTerms(Rb, dq):
    h = np.hstack((Rb.inv().apply([0, 0, -mb * g]), np.zeros(3)))
    B = np.eye(6)
    # FIXME: don't actually need B, since wrenchMap(u) is sort of like B??
	# rot(x) = [cos(x) -sin(x); sin(x) cos(x)]
	# fL = 0.15
	# fR = 0.15
	# # now should be multiplied by u = [ΦL^2; ΦR^2]
	# S = hcat(aeroWrenchAffine(fL), Diagonal([1,-1]) * aeroWrenchAffine(fR)) # 2x2 matrix
	# B = [rot(y[3])[:,2] zeros(2,1); 0 1] * S
    return M, h, B
