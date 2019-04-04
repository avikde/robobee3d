# Import packages.
import cvxpy as cp
import numpy as np

# Generate a random non-trivial quadratic program.
m = 15
n = 10
p = 5
np.random.seed(1)
P = np.random.randn(n, n)
P = P.T@P
q = np.random.randn(n)
G = np.random.randn(m, n)
h = G@np.random.randn(n)
A = np.random.randn(p, n)
b = np.random.randn(p)

# Define and solve the CVXPY problem.
x = cp.Variable(n)
prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(x, P) + q.T@x),
                 [G@x <= h,
                  A@x == b])
# Print result.
prob.solve(solver=cp.MOSEK)
print("\nThe optimal value is", prob.value)

prob.solve(solver=cp.OSQP)
print("\nThe optimal value is", prob.value)

prob.solve(solver=cp.SCS)
print("\nThe optimal value is", prob.value)

prob.solve(solver=cp.ECOS)
print("\nThe optimal value is", prob.value)

prob.solve(solver=cp.ECOS_BB)
print("\nThe optimal value is", prob.value)


