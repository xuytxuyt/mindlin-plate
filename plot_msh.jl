
using ApproxOperator, GLMakie

import Gmsh: gmsh

ndiv = 8
gmsh.initialize()
# gmsh.open("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh")
gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_tri6_"*string(ndiv)*".msh")
# gmsh.open("./msh/circular/circular_quad_"*string(ndiv)*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
elements["Γᵇ"] = getElements(nodes,entities["Γᵇ"])
elements["Γˡ"] = getElements(nodes,entities["Γˡ"])
elements["Γʳ"] = getElements(nodes,entities["Γʳ"])
elements["∂Ω"] = elements["Γᵇ"]∪elements["Γᵗ"]∪elements["Γˡ"]∪elements["Γʳ"]

# elements["Γᵉ"] = getElements(nodes,entities["Γᵉ"])
# elements["∂Ω"] = elements["Γᵇ"]∪elements["Γˡ"]∪elements["Γᵉ"]

nodes_mf_Ω = 𝑿ᵢ[]
for elm in elements["Ω"]
    push!(nodes_mf_Ω,elm.𝓒...)
end
unique!(nodes_mf_Ω)

nodes_mf_Γ = 𝑿ᵢ[]
for elm in elements["∂Ω"]
    push!(nodes_mf_Γ,elm.𝓒...)
end
unique!(nodes_mf_Γ)

# gmsh.finalize()
# gmsh.open("./msh/SquarePlate/SquarePlate_q_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_q_"*string(ndiv)*".msh")
# nodes_mf = get𝑿ᵢ()
# xp = nodes_mf.x
# yp = nodes_mf.y
# zp = nodes_mf.z

f = Figure(backgroundcolor = :transparent)
# axis
ax = Axis3(f[1, 1], perspectiveness = 0.8, aspect = (1,1,0.2), azimuth = 1.5π, elevation = 0.1*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.),)
hidespines!(ax)
hidedecorations!(ax)
# xp = [node.x for node in nodes_mf_Ω] 
# yp = [node.y for node in nodes_mf_Ω] 

# scatter!(ax, xp, yp, marker = :circle, markersize = 8, color = :dodgerblue)

# x2 = [node.x for node in nodes_mf_Γ]
# y2 = [node.y for node in nodes_mf_Γ]
# scatter!(ax, x2, y2, marker = :circle, markersize = 8, color = :black)

# x =  nodes.x
# y = nodes.y
# z = 0
# ps = Point3f.(x,y,z)
# scatter!(ps, 
#     marker=:circle,
#     markersize = 5,
#     color = :black
# )

# elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]

    lines!(x,y,linestyle = :dash, linewidth = 0.5, color = :black)
end

# # boundaries
for elm in elements["∂Ω"]
    ξ¹ = [x.x for x in elm.𝓒]
    ξ² = [x.y for x in elm.𝓒]
    x =  [x.x for x in elm.𝓒]
    y =  [x.y for x in elm.𝓒]
    lines!(x,y,linewidth = 1.5, color = :black)
end
trim!(f.layout)
# save("./png/Circular_"*string(ndiv)*"_msh.png",f)
# save("./png/plate_with_hole_tri3_"*string(ndiv)*"_msh.png",f)
# save("./png/SquarePlate_"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/MorleysAcuteSkewPlate_"*string(ndiv)*"_msh.png",f)
f