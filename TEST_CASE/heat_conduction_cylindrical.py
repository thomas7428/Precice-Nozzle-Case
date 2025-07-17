from nutils import mesh, solver, export
from nutils.expression_v2 import Namespace
import numpy as np


# made up functions for conductivity, Dirichlet BC and Neumann BC (heat flux)
def k(t):
    return 20 + 10 * t


def boundary(coord):
    return 2 * coord**2


def flux(coord):
    return 20 * coord


btype = "std"
degree = 2

# geometry
length = 1.0
inner_diameter = 0.2
outer_diameter = 0.7

# grid generation
grid = [np.linspace(a, b, round((b - a) / size) + 1) for (a, b, size) in [(0, length, 0.01), (inner_diameter, outer_diameter, 0.01)]]
domain, geom = mesh.rectilinear(grid)

# definition of the namespace
ns = Namespace()
ns.x = geom
ns.define_for('x', gradient='d', normal='n', jacobians=('dV', 'dS'))
ns.basis = domain.basis(btype, degree=degree)
ns.t = domain.field('t', btype=btype, degree=degree)
ns.v = domain.field('v', btype=btype, degree=degree)
ns.k = k
ns.boundary = boundary
ns.flux = flux

# variational formulation of non-linear heat equation for 2D cylindrical coordinates. non-linearity results from temperature dependend conductivity
# x_0 --> axial direction
# x_1 --> radial direction
res = domain.integral('(d_1(v / x_1) x_1 d_1(t) k(t) + d_0(v) d_0(t) k(t)) dV' @ ns, degree=degree * 2)

# Neumann boundary condition on bottom boundary calling the flux function
res += domain.boundary['bottom'].integral('v flux(x_0) dS' @ ns, degree=degree * 2)

# Dirichlet boundary condition for top boundary calling the boundary function
sqr = domain.boundary['top'].integral('(t - boundary(x_0))^2 dS' @ ns, degree=degree * 2)

# solving for the constsraints
cons = solver.System(sqr, trial='t').solve_constraints(droptol=1e-15)

# solving the problem
lhs = solver.System(res, trial='t', test='v').solve(constrain=cons, tol=1e-10)

# sampling for output
bezier = domain.sample('bezier', 9)
x, t = bezier.eval(['x_i', 't'] @ ns, lhs)
export.triplot('results.png', x, t, tri=bezier.tri, hull=bezier.hull)
export.vtk('results', bezier.tri, x, Temperature=t)
