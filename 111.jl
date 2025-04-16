
using ApproxOperator
using CairoMakie
using GLMakie

import Gmsh: gmsh

gmsh.initialize()

# boundary
gmsh.open("msh/cloud_0_1.msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

elements["Γ"] = getElements(nodes, entities["Γ"])

# meshfree
gmsh.open("msh/cloud_3.msh")
nodes_mf = get𝑿ᵢ()
x = nodes_mf.x
y = nodes_mf.y
z = nodes_mf.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
s_ = 4
s = s_ .*ones(length(nodes_mf))
push!(nodes_mf,:s=>s)
type = ReproducingKernel{:Linear2D,:○,:CubicSpline}
id = 17
xₑ = nodes_mf[id].x
yₑ = nodes_mf[id].y
inte = 500
xs = LinRange(2.5,9.5,inte)
ys = LinRange(1.5,8,inte)
zs = zeros(inte,inte)
for (i,x) in enumerate(xs)
    for (j,y) in enumerate(ys)
        # r = ((x-xₑ)^2+(y-yₑ)^2)^0.5
        # if r ≤ s_
            indices = sp(x,y,0.0)
            𝓒 = [nodes_mf[i] for i in indices]
            nₛ = length(indices)
            data = Dict([:x=>(2,[x]),:y=>(2,[y]),:z=>(2,[0.]),:𝝭=>(4,zeros(nₛ)),:𝗠=>(0,zeros(6))])
            xₛ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)
            𝓖 = [xₛ]
            a = type(𝓒,𝓖)
            set𝝭!(a)
            N = xₛ[:𝝭]
            zs[i,j] = N[indexin(id,collect(indices))...]
        # else
        #     zs[i,j] = NaN
        # end
    end
end
entities = getPhysicalGroups()
nodes_mf = get𝑿ᵢ()
elements["Γₚ"] = getElements(nodes_mf, entities["Γ"])
nodes_mf_Γ = 𝑿ᵢ[]
for elm in elements["Γₚ"]
    push!(nodes_mf_Γ,elm.𝓒...)
end
unique!(nodes_mf_Γ)


# finite element
gmsh.open("msh/cloud_2.msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements["Ω"] = getElements(nodes, entities["Ω"])

nₚ = length(nodes)
nₑ = length(elements["Ω"])
d = [NaN for i in 1:nₚ]
d[35] = 1.0
d[[29,33,36,37,45]] .= 0.0
vertices = zeros(nₚ,3)
for (i,node) in enumerate(nodes)
    vertices[i,1] = node.x
    vertices[i,2] = node.y
    vertices[i,3] = d[i]
end
faces = zeros(Int,nₑ,3)
for (i,elm) in enumerate(elements["Ω"])
    faces[i,:] .= [x.𝐼 for x in elm.𝓒]
end

# figure
# fig = Figure(backgroundcolor = :transparent) # 1 2
# ax = Axis( # 1 2
fig = Figure(backgroundcolor = :transparent, resolution = (500, 500)) # 3 4
ax = Axis3( # 3 4
    fig[1,1],
    xlabel = " ",
    ylabel = " ",
    xticksvisible = false,
    yticksvisible = false,
    # zticksvisible = false,
    xticklabelsvisible=false,
    yticklabelsvisible=false,
    xgridvisible=false,
    ygridvisible=false,
    # zgridvisible=false,
    # perspectiveness = 0.2,
    backgroundcolor = :transparent,
    
    # aspect = :data,
    aspect = (1,1,0.6), # 0.3(3) 0.2(4)
    azimuth = 1.6π, # 3
    elevation = 0.1π, # 3
    # aspect = AxisAspect(2),
)
hidespines!(ax)
hidedecorations!(ax)

# for elm in elements["Ω"]
#     𝓒 = elm.𝓒
#     x = [node.x for node in 𝓒[[1:end...,1]]]
#     y = [node.y for node in 𝓒[[1:end...,1]]]
    # lines!(x,y, color=:black, linewidth = 1.5) # 1 3
#     # lines!(ax,x,y, color=:lightgrey, linewidth = 1.5) # 2 5
    # lines!(ax,x,y, color=:lightgrey, linewidth = 1.5) # 4
# end
# for elm in elements["Γ"]
#     𝓒 = elm.𝓒
#     x = [node.x for node in 𝓒[1:end]]
#     y = [node.y for node in 𝓒[1:end]]
#     lines!(x,y, color=(:skyblue4,0.0), linewidth = 6) # linewidth = 5(1,2) 3(3,4)
# end


# save("./png/cloud_1.png",fig, px_per_unit = 10.0) # 1


# xₚ = [node.x for node in nodes_mf] # 
# yₚ = [node.y for node in nodes_mf] #
# xₚ = [node.x for node in nodes_mf_Γ] # 2 4
# yₚ = [node.y for node in nodes_mf_Γ] # 2 4
# scatter!(ax, xₚ, yₚ, marker = :circle, markersize = 25, color = :dodgerblue) # 2 4

# x1 = [node.x for node in nodes] # 
# y1 = [node.y for node in nodes] #
# scatter!(ax, x1, y1, marker = :circle, markersize = 20, color = :black) # 2 4

# scatter!(ax, xₚ[id], yₚ[id], marker = :circle, markersize = 20, color = :red)
# p_big = decompose(Point2f, Circle(Point2f(0),1))
# p_small = decompose(Point2f, Circle(Point2f(0),0.99))
# scatter!(ax, xₚ[id], yₚ[id], marker = Polygon(p_big,[p_small]), markersize = 145, color = :dodgerblue) # 2
# scatter!(ax, xₚ, yₚ, marker = Polygon(p_big,[p_small]), markersize = 60, color = :dodgerblue)
# arc!(Point2f(xₚ[id],yₚ[id]),4,0,2π,linewidth = 1.5,color = :salmon)
# for (x,y) in zip(xₚ,yₚ)
#     arc!(Point2f(x,y),1.5,0,2π,linewidth = 1,color = :dodgerblue)
# end
# xₚ = [5.0+(rand()-0.5) for i in 1:8]
# yₚ = [5.0+(rand()-0.5) for i in 1:8]
# scatter!(ax, xₚ, yₚ, marker = :circle, markersize = 20, color = :black)
# for (x,y) in zip(xₚ,yₚ)
#     arc!(Point2f(x,y),1.5,0,2π,linewidth = 1,color = :dodgerblue)
# end
# save("./png/cloud_2.png",fig, px_per_unit = 10.0)
# save("./png/cloud_5.png",fig, px_per_unit = 10.0)
# save("./png/cloud_6.png",fig, px_per_unit = 10.0)
# save("./png/cloud_7.png",fig, px_per_unit = 10.0)

# mesh!(vertices,faces,
#     color=d,
#     # colormap=:jet,
#     colormap=:Spectral,
#     colorrange = (0,1),
#     transparency=false,
#     shading = false,
# )
# save("./png/cloud_3.png",fig, px_per_unit = 20.0) # 3

surface!(xs, ys, zs,
    colormap = :Spectral,
    colorrange = (0.0,1),
    transparency = false,
    shading = false,
)

save("./png/cloud_4.png",fig, px_per_unit = 10.0) # 4

fig
