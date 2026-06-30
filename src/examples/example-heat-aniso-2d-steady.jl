#=
# 2d Advection-Reaction-Diffusion, Dirichlet bc
=#

### If the Finch package has already been added, use this line #########
using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
# if !@isdefined(Finch)
#     include("../Finch.jl");
#     using .Finch
# end
##########################################################################

initFinch("heat-anisotropic-steady-2d.jl");
useLog("heat-anisotropic-steady-2d.log", level=3)

n = 32;
ord = 1;

domain(2)

functionSpace(order=ord)

mesh(TRIMESH, elsperdim=n)

u = variable("u")
testSymbol("v")

boundary(u, 1, DIRICHLET, 0)


# Material properties - fiber-reinforced composite
k_parallel = 100.0    # Conductivity along fibers
k_perpendicular = 1.0 # Conductivity perpendicular to fibers
theta_deg = 60.0      # Fiber orientation angle

# Heat source
Q0 = 1.0

# -------------------------------------------------------------------
# Compute rotated conductivity tensor
# -------------------------------------------------------------------
theta = deg2rad(theta_deg)
c, s = cos(theta), sin(theta)
R = [c -s; s c]

D_principal = [k_parallel 0.0; 0.0 k_perpendicular]

# Transform to global coordinates: K = R * D * R^T
K = R * D_principal * R'

# Write the weak form 
coefficient("Q", -Q0)
coefficient("Dtensor", K, type=TENSOR)


weakForm(u, "dot(matvec(Dtensor, grad(u)), grad(v)) + Q*v")

exportCode("heat-anisotropic-steady-2dcode");
# importCode("heat-anisotropic-steady-2dcode");

solve(u);

outputValues(u, "heat-anisotropic-steady-2doutput"; format="vtk", ascii=false)
finalizeFinch()
