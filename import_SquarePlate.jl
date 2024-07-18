
using Tensors, BenchmarkExample, Statistics, DelimitedFiles
import Gmsh: gmsh
function import_SquarePlate(filename::String)
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
    # elements["𝐴"] = getElements(nodes, entities["𝐴"], integrationOrder)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    # gmsh.finalize()
    return elements, nodes
end
function import_SquarePlate_p(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    integrationOrder = 12
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z

    type = PiecewisePolynomial{:6}
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getPiecewiseElements(entities["Ω"],type, integrationOrder)
    # gmsh.finalize()
    return elements, nodes
end

function import_SquarePlate_mix(filename1::String,filename2::String)
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

    gmsh.open(filename2)
    nodes_s = get𝑿ᵢ()
    xˢ = nodes_s.x
    yˢ = nodes_s.y
    zˢ = nodes_s.z
    # s = 2.5/ndivs*ones(length(nodes_s))
    Ω = getElements(nodes_s, entities["Ω"])
    s, var𝐴 = cal_area_support(Ω)
    s = 2.5*s*ones(length(nodes_s))
    push!(nodes_s,:s₁=>s,:s₂=>s,:s₃=>s)
    type = ReproducingKernel{:Linear2D,:□,:CubicSpline}
    sp = RegularGrid(xˢ,yˢ,zˢ,n = 1,γ = 2)

    gmsh.open(filename1)
    elements["Ωˢ"] = getElements(nodes_s, entities["Ω"], type, integrationOrder, sp)
    nₘ=21
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    push!(elements["Ωˢ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωˢ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    # gmsh.finalize()
    return elements, nodes, nodes_s
end
function import_SquarePlate_quad_RI(filename1::String,filename2::String)
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
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], integrationOrder,normal=true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], integrationOrder,normal=true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], integrationOrder,normal=true)
    elements["Γʳ"] = getElements(nodes, entities["Γʳ"], integrationOrder,normal=true)
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
prescribeForSSNonUniformLoading = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵇ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γᵗ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γˡ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γʳ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γᵇ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γᵗ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γʳ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γᵗ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γˡ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γʳ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Ω"],:q=>(x,y,z)->F(x,y))
end

prescribeForSSUniformLoading = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵇ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γˡ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γʳ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γᵇ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γʳ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Γˡ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Γʳ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Ω"],:q=>(x,y,z)->F)
end

prescribeForSimpleSupported = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵇ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γᵗ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γˡ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γʳ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γᵇ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γᵗ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γʳ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γᵗ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γˡ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γʳ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Ω"],:q=>(x,y,z)->F(x,y))
    prescribe!(elements["Γᵇ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γᵗ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γˡ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γʳ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γᵇ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γᵗ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γˡ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γʳ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γᵇ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γᵗ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γˡ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γʳ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
end
prescribeForSSUniformLoading = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵇ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γˡ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γʳ"],:g=>(x,y,z)->0.0)
    prescribe!(elements["Γᵇ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γʳ"],:θ₁=>(x,y,z)->0.0)
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Γᵗ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Γˡ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Γʳ"],:θ₂=>(x,y,z)->0.0)
    prescribe!(elements["Ω"],:q=>(x,y,z)->F)
end

prescribeForCantilever = quote
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωᵍ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠)

    prescribe!(elements["Γᵇ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γᵗ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γˡ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γʳ"],:g=>(x,y,z)->w(x,y))
    prescribe!(elements["Γᵇ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γᵗ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γʳ"],:θ₁=>(x,y,z)->θ₁(x,y))
    prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γᵗ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γˡ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Γʳ"],:θ₂=>(x,y,z)->θ₂(x,y))
    prescribe!(elements["Ω"],:q=>(x,y,z)->F(x,y))
    prescribe!(elements["Γᵇ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γᵗ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γˡ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γʳ"],:M₁₁=>(x,y,z)->M₁₁(x,y))
    prescribe!(elements["Γᵇ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γᵗ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γˡ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γʳ"],:M₁₂=>(x,y,z)->M₁₂(x,y))
    prescribe!(elements["Γᵇ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γᵗ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γˡ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γʳ"],:M₂₂=>(x,y,z)->M₂₂(x,y))
    prescribe!(elements["Γᵇ"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γᵗ"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γˡ"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γʳ"],:Q₁=>(x,y,z)->Q₁(x,y))
    prescribe!(elements["Γᵇ"],:Q₂=>(x,y,z)->Q₂(x,y))
    prescribe!(elements["Γᵗ"],:Q₂=>(x,y,z)->Q₂(x,y))
    prescribe!(elements["Γˡ"],:Q₂=>(x,y,z)->Q₂(x,y))
    prescribe!(elements["Γʳ"],:Q₂=>(x,y,z)->Q₂(x,y))
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