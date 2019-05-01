import fenics as fe
import dolfin
import numpy as np
from dolfin.fem.norms import errornorm
from dolfin.common.plotting import plot
import matplotlib.pyplot as plt
import sys

EPSILON = 1.0e-14
DEG = 2

mesh = fe.Mesh('step.xml')

# Control pannel
MODEL = False # flag to use SA model
b = fe.Expression(('0', '0'), degree=DEG) # forcing
nu = fe.Constant(.05)
rho = fe.Constant(1)
RE = 50
lmx = .3 # mixing length 
dt = 0.01
# Re = 10 / 1e-4 = 1e5

V   = fe.VectorElement("Lagrange", mesh.ufl_cell(), 2)
P   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
NU   = fe.FiniteElement("Lagrange", mesh.ufl_cell(), 1)
if MODEL: M   = fe.MixedElement([V, P, NU])
else: M   = fe.MixedElement([V, P])
W   = fe.FunctionSpace(mesh, M)

W0 = fe.Function(W)
We = fe.Function(W)
u0, p0 = fe.split(We)
#u0 = fe.Function((W0[0], W0[1]), 'Velocity000023.vtu')
#p0 = fe.Function(W0[2])
v, q = fe.TestFunctions(W)
#u, p = fe.split(W0)
u,p = (fe.as_vector((W0[0], W0[1])), W0[2])




#-------------------------------------------------------
# Defining essential/Dirichlet boundary conditions
# Step 1: Identify all boundary segments forming Gamma_d
#-------------------------------------------------------
#   (-3., 2.5)
#  |
#  |
#  |_______(0, 1.)
#    bc1  |
#      bc2|__________(3., ,0.)
#        (0,0) bc3

# surface before step
def dbc1(x, on_boundary):
    return on_boundary and np.abs(x[1] - 1.) < EPSILON and x[0] < EPSILON

# surface on step side
def dbc2(x, on_boundary):
    return on_boundary and x[0] < EPSILON and x[1] < 1.0

# surface after step
def dbc3(x, on_boundary):
    return on_boundary and x[1] < EPSILON and x[0] > - 1 * EPSILON

# inflow
def dbc_inflow(x, on_boundary):
    return on_boundary and np.abs(x[0] + 0.2) < EPSILON

# outlet
def dbc_outflow(x, on_boundary):
    return on_boundary and np.abs(x[0] - 10) < EPSILON

# top
def dbc_top(x, on_boundary):
    return on_boundary and np.abs(x[1] - 2.5) < EPSILON

#--------------------------------------------------------
# Defining essential/Dirichlet boundary conditions
# Step 2: Defining what the boundary values will be (u_D)
#--------------------------------------------------------
uD_X0 = fe.Expression(('0', '0'), degree=DEG)
uD_X1 = fe.Expression(('0', '0'), degree=DEG)
uD_Y0 = fe.Expression(('0', '0'), degree=DEG)
uD_Y1 = fe.Expression(('%s * pow((x[1] - 1), 2)' % RE, '0'), degree=DEG)
#uD_Y1 = fe.Expression(('%s' % RE, '0'), degree=DEG)
bc_p  = fe.Constant(('0'))
bc_v  = fe.Constant(('0'))
bc_vf  = fe.Constant(3 * nu)

bc_1 = fe.DirichletBC(W.sub(0), uD_X0, dbc1)
bc_2 = fe.DirichletBC(W.sub(0), uD_Y0, dbc2)
bc_3 = fe.DirichletBC(W.sub(0), uD_X1, dbc3)
#bc_inflow = fe.DirichletBC(W.sub(0), fe.Expression(('%s * pow((x[1] - 1) / 2.5, 2) + 0.01 * sin(20 * pi * x[1])' % RE, '0'), degree=DEG), dbc_inflow)
bc_inflow = fe.DirichletBC(W.sub(0), fe.Expression(('%s * pow((x[1] - 1) / 2.5, 2)' % RE, '0'), degree=DEG), dbc_inflow)
#bc_inflow = fe.DirichletBC(W.sub(0), fe.Constant((RE, '0')), dbc_inflow)
bc_p = fe.DirichletBC(W.sub(1), bc_p, dbc_top)

def Max(a, b): return (a + b + abs(a-b)) / 2.
def Min(a, b): return (a + b - abs(a-b)) / 2.


ns_conv = fe.inner(v, fe.grad(u)*u)*fe.dx
ns_press = p * fe.div(v) * fe.dx
#s = fe.grad(u) + fe.grad(u).T
sij = 0.5 * (fe.grad(u) + fe.grad(u).T)
S = fe.sqrt(EPSILON + 2.*fe.inner(0.5*(fe.grad(u)+fe.grad(u).T), 0.5*(fe.grad(u)+fe.grad(u).T)))
lmx = 0.01
nu_tv = lmx ** 2. * S

#nu_tv = 0.5 * fe.inner(sij, sij) ** 0.5
#nu_tv = lmx * (2 * fe.inner(sij, sij)) ** (0.5)
#ns_tv = fe.inner((nu_tv) * fe.grad(v), fe.grad(u)) * fe.dx
ns_visc = (nu + nu_tv) * fe.inner(fe.grad(v), fe.grad(u)) * fe.dx
#ns_visc = nu * fe.inner(fe.grad(v), fe.grad(u)) * fe.dx
ns_conti = q * fe.div(u) * fe.dx
ns_forcing = fe.dot(v, b)*fe.dx

NS = ns_conv + ns_press + ns_visc - ns_conti + ns_forcing
#NS = ns_conv + ns_press + ns_tv + ns_visc + ns_conti + ns_forcing

N = 5
#fe.parameters["form_compiler"]["quadrature_degree"] = N

weakForm  = (1.0 / dt) * fe.inner(u-u0, v) * fe.dx + NS

STAB = False
he  = fe.CellDiameter(mesh)
tau = (1.0/3.0)*(he*he)/(4.0 * nu * rho)

res = - tau * fe.inner(fe.dot(u, fe.grad(v)), fe.grad(u)*u) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)),  fe.grad(p)) * fe.dx
#res += - tau * fe.inner(fe.dot(u, fe.grad(v)),  -1 * fv1 * nu_trial * fe.div(fe.grad(u))) * fe.dx #TODO - update residual term for new turbulence model
res += - tau * fe.inner(fe.dot(u, fe.grad(v)), -1 * nu * fe.div(fe.grad(u))) * fe.dx 
res += - tau * fe.inner(fe.dot(u, fe.grad(q)), fe.div(u)) * fe.dx
res += - tau * fe.inner(fe.dot(u, fe.grad(v)), -1 * b) * fe.dx 

if STAB: 
   weakForm += res
stab     = -tau*fe.inner(fe.grad(q), fe.grad(p))*fe.dx
#weakForm = weakForm + stab

dW  = fe.TrialFunction(W)
dFdW = fe.derivative(weakForm, W0, dW)

if MODEL: bcSet   = [bc_1, bc_2, bc_3, bc_inflow, bc_p, bc_v_x0, bc_v_x1, bc_v_y1, bc_v_in, bc_v_top]
else: bcSet   = [bc_1, bc_2, bc_3, bc_inflow, bc_p]
problem = fe.NonlinearVariationalProblem(weakForm, W0, bcSet, J=dFdW)

solver = fe.NonlinearVariationalSolver(problem)

prm = solver.parameters


t = 0.0
t_end = 5
pFile = fe.File('Pressure.pvd')
uFile = fe.File('Velocity.pvd')
vFile = fe.File('Vorticity.pvd')
wFile = fe.File('W.pvd')
T = fe.FunctionSpace(mesh, 'CG', 1)
#solver.solve()
ii = 0
while t < t_end:
    print("t =",t)
    solver.solve()
    u1,p1 = W0.split()
    We.assign(W0)
    if True:
      uFile << u1
      pFile << p1
    #wFile << W0
      #omega = fe.curl(u)
      #vFile << fe.project(fe.inner(omega, omega), T)
    t += dt
    ii += 1

u, p = W0.split()

#-------------------------------------------------
# Save this solution to a file for post-processing
#-------------------------------------------------
T = fe.TensorFunctionSpace(mesh, 'CG', 1)
wFile = fe.File('Sij.pvd')
wFile << fe.project(sij, T)
T = fe.FunctionSpace(mesh, 'CG', 1)
wFile = fe.File('S.pvd')
wFile << fe.project(lmx**2.*S, T)
T = fe.VectorFunctionSpace(mesh, 'CG', 1)
vtkFile = fe.File('CFL.pvd')
vtkFile << fe.project(u / he, T)

T = fe.TensorFunctionSpace(mesh, 'CG', 1)
vtkFile = fe.File('tau.pvd')
vtkFile << fe.project(nu * (fe.grad(u) + fe.grad(u).T), T)

T = fe.FunctionSpace(mesh, 'CG', 1)
wFile = fe.File('S.pvd')
wFile << fe.project(lmx**2.*S, T)
