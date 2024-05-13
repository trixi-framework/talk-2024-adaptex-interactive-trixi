using HOHQMesh

# plotting backend
using GLMakie

project = newProject("Mountain", "output_mesh")  # project name, output folder

mountainHeight = 400      # height of montain 400 m
mountainHalfWidth = 1000  # half width of mountain 1 km
x_ll = -12000.0           # domain width 40 km
x_lr =  28000.0
width = x_lr - x_ll

# lower left and right y-values according to parametrization
y_ll = mountainHeight / ( 1 + (x_ll / mountainHalfWidth)^2 )
y_lr = mountainHeight / ( 1 + (x_lr / mountainHalfWidth)^2 )

# parametrized bottom curve, t ∈ [0,1]
xEqn = "x(t) = $x_ll + t * $width"
yEqn = "y(t) = $mountainHeight / ( 1 + (($x_ll + t * $width) / $mountainHalfWidth)^2 )"
zEqn = "z(t) = 0.0"
bottom = newParametricEquationCurve("bottom", xEqn, yEqn, zEqn)

# spline curve above
spline_data = [ [0.0   x_lr                 y_lr          0.0]
                [0.2   x_ll + 0.9 * width   0.4 * width   0.0]
                [0.8   x_ll + 0.1 * width   0.4 * width   0.0]
                [1.0   x_ll                 y_ll          0.0] ]

dome = newSplineCurve("dome", size(spline_data)[1], spline_data)

# use curves as outer boundary
for curve in [dome, bottom]
    addCurveToOuterBoundary!(project, curve)
end

# set the desired resolution in x and y direction
addBackgroundGrid!(project, [2000.0, 2000.0, 0.0])

# have a look
plotProject!(project, MODEL+GRID)

# set polynomial degree for curved boundaries
setPolynomialOrder!(project, 3)

# change output format to ABAQUS
setMeshFileFormat!(project, "ABAQUS")

# generate
@info "Press enter to generate the mesh and update the plot."
readline()
generate_mesh(project)

@info "Press enter to quit."
readline()
