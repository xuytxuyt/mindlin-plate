
using ApproxOperator, GLMakie

import Gmsh: gmsh

ndiv =  32
gmsh.initialize()
# gmsh.open("./msh/cook_tri3_"*string(ndiv)*".msh")
# gmsh.open("./msh/plate_with_hole_tri3_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad8_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_tri6_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_quad_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_quad8_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_tri6_"*string(ndiv)*".msh")
# gmsh.open("./msh/circular/circular_"*string(ndiv)*".msh")
gmsh.open("./msh/circular/circular_quad_"*string(ndiv)*".msh")
# gmsh.open("./msh/circular/circular_quad8_"*string(ndiv)*".msh")
# gmsh.open("./msh/circular/circular_tri6_"*string(ndiv)*".msh")
# gmsh.open("./msh/cantilever_"*string(ndiv)*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()

elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
# elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
elements["Γᵇ"] = getElements(nodes,entities["Γᵇ"])
elements["Γˡ"] = getElements(nodes,entities["Γˡ"])
# elements["Γʳ"] = getElements(nodes,entities["Γʳ"])
# elements["∂Ω"] = elements["Γᵇ"]∪elements["Γᵗ"]∪elements["Γˡ"]∪elements["Γʳ"]

elements["Γᵉ"] = getElements(nodes,entities["Γᵉ"])
elements["∂Ω"] = elements["Γᵇ"]∪elements["Γˡ"]∪elements["Γᵉ"]

# elements["Γᵗ"] = getElements(nodes,entities["Γᵗ"])
# elements["Γᵍ"] = getElements(nodes,entities["Γᵍ"])
# elements["∂Ω"] = elements["Γᵍ"]∪elements["Γᵗ"]

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
# gmsh.open("./msh/SquarePlate/SquarePlate_28.msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_q_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_quad_q_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_quad8_q_"*string(ndiv)*".msh")
# gmsh.open("./msh/SquarePlate/SquarePlate_tri6_q_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad_q_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad8_q_"*string(ndiv)*".msh")
# gmsh.open("./msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_tri6_q_"*string(ndiv)*".msh")
# nodes_mf = get𝑿ᵢ()
# xp = nodes_mf.x
# yp = nodes_mf.y
# zp = nodes_mf.z

f = Figure()
# axis
ax = Axis3(f[1, 1], perspectiveness = 0.8, aspect = :data, azimuth = 1.5π, elevation = 0.5*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.),)
hidespines!(ax)
hidedecorations!(ax)
# xp = [node.x for node in nodes_mf] 
# yp = [node.y for node in nodes_mf] 

# scatter!(ax, xp, yp, marker = :circle, markersize = 30, color = :dodgerblue)

# x2 = [node.x for node in nodes_mf_Γ]
# y2 = [node.y for node in nodes_mf_Γ]
# scatter!(ax, x2, y2, marker = :circle, markersize = 30, color = :black)

x =  nodes.x
y = nodes.y
z = 0
ps = Point3f.(x,y,z)
scatter!(ps, 
    marker=:circle,
    markersize = 10,
    color = :black
)

# elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3]]]
    y = [x.y for x in elm.𝓒[[1,2,3]]]

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
# save("./png/cantilever_"*string(ndiv)*"_msh.png",f)
# save("./png/Circular_"*string(ndiv)*"_msh.png",f)
save("./png/Circular_quad4_"*string(ndiv)*"_msh.png",f)
# save("./png/Circular_quad8_"*string(ndiv)*"_msh.png",f)
# save("./png/Circular_tri6_"*string(ndiv)*"_msh.png",f)
# save("./png/plate_with_hole_"*string(ndiv)*"_msh.png",f)
# save("./png/cook_"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/Square_"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/SquarePlate_"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/SquarePlate_tri6_"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/SquarePlate_quad4_"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/SquarePlate_quad8_"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/MorleysAcuteSkewPlate_tri3"*string(ndiv)*"_msh.png",f, px_per_unit = 10.0)
# save("./png/MorleysAcuteSkewPlate_quad4_"*string(ndiv)*"_msh.png",f)
# save("./png/MorleysAcuteSkewPlate_quad8_"*string(ndiv)*"_msh.png",f)
# save("./png/MorleysAcuteSkewPlate_tri6_"*string(ndiv)*"_msh.png",f)
f