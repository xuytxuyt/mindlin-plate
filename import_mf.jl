
using Tensors, BenchmarkExample, Statistics, DelimitedFiles
import Gmsh: gmsh
function import_mf(filename::String)
    gmsh.initialize()
    gmsh.open(filename)
    integrationOrder = 2     # Tri3
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    xˢ = nodes.x
    yˢ = nodes.y
    zˢ = nodes.z
    s = 1.5/ndiv*ones(length(nodes))
    push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xˢ,yˢ,zˢ,n = 1,γ = 2)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], type, integrationOrder, sp)
    nₘ=21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # gmsh.finalize()
    return elements, nodes,  sp, type
end