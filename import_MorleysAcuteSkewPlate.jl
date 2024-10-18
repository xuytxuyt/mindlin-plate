
using Tensors, BenchmarkExample
import Gmsh: gmsh
function import_MorleysAcuteSkewPlate(filename::String)
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
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], integrationOrder,normal=true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], integrationOrder,normal=true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], integrationOrder,normal=true)
    elements["Γʳ"] = getElements(nodes, entities["Γʳ"], integrationOrder,normal=true)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], integrationOrder)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes
end
function import_MorleysAcuteSkewPlate_mix(filename1::String,filename2::String)
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
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], integrationOrder,normal=true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], integrationOrder,normal=true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], integrationOrder,normal=true)
    elements["Γʳ"] = getElements(nodes, entities["Γʳ"], integrationOrder,normal=true)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], integrationOrder)

    gmsh.open(filename2)
    entities = getPhysicalGroups()
    nodes_s = get𝑿ᵢ()
    xˢ = nodes_s.x
    yˢ = nodes_s.y
    zˢ = nodes_s.z
    s = 1.5*100/ndivs*ones(length(nodes_s))
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
    return elements, nodes, nodes_s, Ω
end

prescribeForSSUniformLoading = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵇ"],:g=>(x,y,z)->w)
    prescribe!(elements["Γᵗ"],:g=>(x,y,z)->w)
    prescribe!(elements["Γˡ"],:g=>(x,y,z)->w)
    prescribe!(elements["Γʳ"],:g=>(x,y,z)->w)
    prescribe!(elements["Γᵇ"],:θ₁=>(x,y,z)->θ₁)
    prescribe!(elements["Γᵗ"],:θ₁=>(x,y,z)->θ₁)
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->θ₁)
    prescribe!(elements["Γʳ"],:θ₁=>(x,y,z)->θ₁)
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->θ₂)
    prescribe!(elements["Γᵗ"],:θ₂=>(x,y,z)->θ₂)
    prescribe!(elements["Γˡ"],:θ₂=>(x,y,z)->θ₂)
    prescribe!(elements["Γʳ"],:θ₂=>(x,y,z)->θ₂)
    prescribe!(elements["Ω"],:q=>(x,y,z)->F)

end