#=
# 2D Poisson, Dirichlet bc
# CG, Linear element, triangle mesh
# Linear system solved with PETSc using a direct LU factorization (PCLU).
# Serial only: no partitioning or MPI.
=#

### If the Finch package has already been added, use this line #########
# using Finch # Note: to add the package, first do: ]add "https://github.com/paralab/Finch.git"

### If not, use these four lines (working from the examples directory) ###
if !@isdefined(Finch)
    include("../Finch.jl");
    using .Finch
end
##########################################################################

initFinch("poisson2dpetsc");

useLog("poisson2dpetsclog", level=3)

# Solve with PETSc: KSP type preonly + PC type lu
linAlgOptions(pc="LU")
usePetsc()

# Set up the configuration
domain(2) # dimension
functionSpace(order=1) # linear elements

mesh(TRIMESH, elsperdim=32)   # uniform triangle mesh of the unit square

u = variable("u")              # make a scalar variable with symbol u
testSymbol("v")                # sets the symbol for a test function

boundary(u, 1, DIRICHLET, 0)  # boundary condition for BID 1 is Dirichlet with value 0

# Write the weak form for -laplacian(u) = f
coefficient("f", "2*pi*pi*sin(pi*x)*sin(pi*y)")
weakForm(u, "dot(grad(u), grad(v)) - f*v")

solve(u);

outputValues(u, "poisson2dpetscoutput"; format="vtk", ascii=false)

finalizeFinch()

# exact solution is sin(pi*x)*sin(pi*y)
# check error
allerr = zeros(size(Finch.finch_state.grid_data.allnodes,2));

for i=1:size(Finch.finch_state.grid_data.allnodes,2)
    x = Finch.finch_state.grid_data.allnodes[1,i];
    y = Finch.finch_state.grid_data.allnodes[2,i];
    exact = sin(pi*x)*sin(pi*y);
    allerr[i] = abs(u.values[i] - exact);
end
maxerr = maximum(abs, allerr);
println("max error = "*string(maxerr));
