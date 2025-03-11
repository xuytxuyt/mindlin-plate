using ApproxOperator, GLMakie
import Gmsh: gmsh

gmsh.initialize()

include("import_mf.jl")

gmsh.open("msh/mf_5.msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
s = 1.5/4 .*ones(length(nodes))
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
# type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
# push!(nodes,:s=>s)
# type = ReproducingKernel{:Linear2D,:○,:CubicSpline}

id = 21
xₑ = nodes[id].x
yₑ = nodes[id].y
inte = 500
xs = LinRange(0.0,1.0,inte)
ys = LinRange(0.0,1.0,inte)
zs = zeros(inte,inte)
for (i,x) in enumerate(xs)
   for (j,y) in enumerate(ys)
      indices = sp(x,y,0.0)
      𝓒 = [nodes[i] for i in indices]
      nₛ = length(indices)
      data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.]),:𝝭=>(4,zeros(nₛ)),:𝗠=>(0,zeros(6))])
      xₛ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)
      𝓖 = [xₛ]
      a = type(𝓒,𝓖)
      set𝝭!(a)
      N = xₛ[:𝝭]
      zs[i,j] = N[indexin(id,collect(indices))...]
   end
end
fig = Figure() 

ax = Axis3( 
    fig[1,1],
    xlabel = " ",
    ylabel = " ",
    zlabel = " ",
   #  zticks = 0:0.4:0.8 ,
    xticksvisible = false,
    yticksvisible = false,
   #  zticksvisible = false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
   #  zticklabelsvisible=false,
    xgridvisible=false,
    ygridvisible=false,
    zgridvisible=false,
    zticklabelsize =30,
    backgroundcolor = :transparent,
   aspect = (1,1,0.8),
)  
# hidespines!(ax)
# hidedecorations!(ax)
limits!(ax, (0, 1), (0, 1), (-0.1, 1))
surface!(xs, ys, zs,
    colormap = :Spectral,
    colorrange = (0,1),
    transparency = false,
    shading = false,
)
save("./png/shapefunction11.png",fig, px_per_unit = 3.0) 
 fig

     
