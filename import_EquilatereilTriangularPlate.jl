
using Tensors, BenchmarkExample, Statistics, DelimitedFiles
import Gmsh: gmsh
function import_EquilatereilTriangularPlate(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationOrder = 2     # Tri3
    # integrationOrder = 3     # Quad4 
    integrationOrder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γ₁"] = getElements(nodes, entities["Γ₁"], integrationOrder,normal=true)
    elements["Γ₂"] = getElements(nodes, entities["Γ₂"], integrationOrder,normal=true)
    elements["Γ₃"] = getElements(nodes, entities["Γ₃"], integrationOrder,normal=true)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes
end


function import_EquilatereilTriangularPlate_mix(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    integrationOrder = 2      # Tri3
    # integrationOrder = 3      # Quad4
    integrationOrder_Ωᵍ = 10
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], integrationOrder)
    elements["Ωᵍ"] = getElements(nodes, entities["Ω"], integrationOrder_Ωᵍ)
    elements["Γ₁"] = getElements(nodes, entities["Γ₁"], integrationOrder,normal=true)
    elements["Γ₂"] = getElements(nodes, entities["Γ₂"], integrationOrder,normal=true)
    elements["Γ₃"] = getElements(nodes, entities["Γ₃"], integrationOrder,normal=true)

    gmsh.open(filename2)
    nodes_s = get𝑿ᵢ()
    xˢ = nodes_s.x
    yˢ = nodes_s.y
    zˢ = nodes_s.z
    s = 1.5*20/3^0.5/ndivs*ones(length(nodes_s))
    Ω = getElements(nodes_s, entities["Ω"])
    push!(nodes_s,:s₁=>s,:s₂=>s,:s₃=>s)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xˢ,yˢ,zˢ,n = 1,γ = 2)

    gmsh.open(filename1)
    elements["Ωᵍˢ"] = getElements(nodes_s, entities["Ω"],type, integrationOrder_Ωᵍ, sp)
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], type, integrationOrder, sp)
    nₘ=21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωˢ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωᵍˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍˢ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # gmsh.finalize()
    return elements, nodes, nodes_s, Ω, sp, type
end

prescribeFor = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ₁"], :𝝭=>:𝑠)
    push!(elements["Γ₂"], :𝝭=>:𝑠)
    push!(elements["Γ₃"], :𝝭=>:𝑠)

    prescribe!(elements["Γ₁"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γ₂"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γ₃"],:g=>(x,y,z)->0.0)
    # prescribe!(elements["Γ₁"],:θ₁=>(x,y,z)->0.0)
    # prescribe!(elements["Γ₂"],:θ₁=>(x,y,z)->0.0)
    # prescribe!(elements["Γ₃"],:θ₁=>(x,y,z)->0.0)
    # prescribe!(elements["Γ₁"],:θ₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γ₂"],:θ₂=>(x,y,z)->0.0)
    # prescribe!(elements["Γ₃"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Ω"],:q=>(x,y,z)->q)
end

function cal_area_support(elms::Vector{ApproxOperator.AbstractElement})
    𝐴s = zeros(length(elms))
    for (i,elm) in enumerate(elms)
        x₁ = elm.𝓒[1].x
        y₁ = elm.𝓒[1].y
        x₂ = elm.𝓒[2].x
        y₂ = elm.𝓒[2].y
        x₃ = elm.𝓒[3].x
        y₃ = elm.𝓒[3].y
        𝐴s[i] = 0.5*(x₁*y₂ + x₂*y₃ + x₃*y₁ - x₂*y₁ - x₃*y₂ - x₁*y₃)
    end
    avg𝐴 = mean(𝐴s)
    var𝐴 = var(𝐴s)
    s = (4/3^0.5*avg𝐴)^0.5
    return s, var𝐴
end