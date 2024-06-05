
using Tensors, BenchmarkExample
import Gmsh: gmsh

function import_patch_test_fem(filename::String)
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
    elements["Γ₄"] = getElements(nodes, entities["Γ₄"], integrationOrder,normal=true)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes
end

function import_patch_test_quad_RI(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    integrationOrder = 3
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
    elements["Γ₄"] = getElements(nodes, entities["Γ₄"], integrationOrder,normal=true)

    gmsh.open(filename2)
    integrationOrder_Ωˢ = 0
    entities = getPhysicalGroups()
    nodes_s = get𝑿ᵢ()
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], integrationOrder_Ωˢ)
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes, nodes_s
end

function import_patch_test_mix(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    # integrationOrder = 2      # Tri3
    integrationOrder = 3      # Quad4
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
    elements["Γ₄"] = getElements(nodes, entities["Γ₄"], integrationOrder,normal=true)

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_s = get𝑿ᵢ()
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], integrationOrder)
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes, nodes_s
end

function import_patch_test_quad_mix(filename1::String,filename2::String)
    gmsh.initialize()
    gmsh.open(filename1)

    integrationOrder = 3
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
    elements["Γ₄"] = getElements(nodes, entities["Γ₄"], integrationOrder,normal=true)

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_s = get𝑿ᵢ()
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], integrationOrder)
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes, nodes_s
end
prescribeForFem = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ₁"], :𝝭=>:𝑠)
    push!(elements["Γ₂"], :𝝭=>:𝑠)
    push!(elements["Γ₃"], :𝝭=>:𝑠)
    push!(elements["Γ₄"], :𝝭=>:𝑠)

    prescribe!(elements["Γ₁"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γ₂"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γ₃"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γ₄"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γ₁"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γ₂"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γ₃"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γ₄"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γ₁"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γ₂"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γ₃"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γ₄"],:θ₂=>(x,y,z)->θ₂(x,y))

    prescribe!(elements["Γ₁"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γ₂"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γ₃"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γ₄"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γ₁"],:Q₂=>(x,y,z)->Q₂(x,y))
    prescribe!(elements["Γ₂"],:Q₂=>(x,y,z)->Q₂(x,y))
    prescribe!(elements["Γ₃"],:Q₂=>(x,y,z)->Q₂(x,y))
    prescribe!(elements["Γ₄"],:Q₂=>(x,y,z)->Q₂(x,y))
    prescribe!(elements["Γ₁"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γ₂"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γ₃"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γ₄"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γ₁"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γ₂"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γ₃"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γ₄"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γ₁"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γ₂"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γ₃"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γ₄"],:M₂₂=>(x,y,z)->M₂₂(x,y))

    # prescribe!(elements["Ω"],:u=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Ω"],:u=>(x,y,z)->w(x,y))
    # prescribe!(elements["Ω"],:∂u∂x=>(x,y,z)->w₁(x,y))
    # prescribe!(elements["Ω"],:∂u∂y=>(x,y,z)->w₂(x,y))
    # prescribe!(elements["Ω"],:M₁ᵢᵢ=>(x,y,z)->M₁₁₁(x,y)+M₁₂₂(x,y))
    # prescribe!(elements["Ω"],:M₂ᵢᵢ=>(x,y,z)->M₁₂₁(x,y)+M₂₂₂(x,y))
    prescribe!(elements["Ω"],:q=>(x,y,z)->-Q₁₁(x,y)-Q₂₂(x,y))
    prescribe!(elements["Ω"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Ω"],:Q₂=>(x,y,z)->Q₂(x,y))

end
