using ApproxOperator, JLD, XLSX, LinearAlgebra

import BenchmarkExample: BenchmarkExample

include("import_EquilatereilTriangularPlate.jl")
ndiv = 10
elements, nodes = import_EquilatereilTriangularPlate("msh/triangle_"*string(ndiv)*".msh");

nₚ = length(nodes)

E = 3e6
ν = 0.3
h = 10
L = 10.0
q = 1.0
D = E*h^3/12/(1-ν^2)
G = 5/6*h*E/2/(1+ν)
# w(x,y) = q/(4*L*D)*(x^3-3*y^2*x-L*(x^2+y^2)+4/27*L^3)*((4/9*L^2-x^2-y^2)/16+D/G)
w(x,y) = q/(64*L*D)*(x^3-3*y^2*x-L*(x^2+y^2)+4/27*L^3)*(4/9*L^2-x^2-y^2)-q/(64*L)*(-14*x^3-64/27*L^3+48*y^2*x+16*L*y^2+16*L*x^2)/G

eval(prescribeFor)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])

ops = [
    Operator{:∫κMγQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫γQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wqdΩ}(),
    Operator{:∫vwdΓ}(:α=>1e13*E),
    Operator{:∫vθ₁dΓ}(:α=>1e13*E),
    Operator{:∫vθ₂dΓ}(:α=>1e13*E),
    Operator{:L₂_ThickPlate_w}(:E=>E,:ν=>ν),
]
k = zeros(3*nₚ,3*nₚ)
kᵇ = zeros(3*nₚ,3*nₚ)
kˢ = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)

# ops[1](elements["Ω"],k)
ops[2](elements["Ω"],kᵇ)
ops[3](elements["Ω"],kˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Γ₁"],k,f)
ops[5](elements["Γ₂"],k,f)
ops[5](elements["Γ₃"],k,f)
# ops[6](elements["Γ₁"],k,f)
# ops[6](elements["Γ₂"],k,f)
# ops[6](elements["Γ₃"],k,f)
# ops[7](elements["Γ₁"],k,f)
# ops[7](elements["Γ₂"],k,f)
# ops[7](elements["Γ₃"],k,f)

d = (kᵇ+kˢ+k)\f
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
set𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωᵍ"])
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->w(x,y))
# prescribe!(elements["Ωᵍ"],:θ₁=>(x,y,z)->θ₁(x,y))
# prescribe!(elements["Ωᵍ"],:θ₂=>(x,y,z)->θ₂(x,y))
L₂ = ops[8](elements["Ωᵍ"])
a = log10(L₂)
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/SquarePlate.xlsx", mode="rw") do xf
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["B"*string(ind)] = log10(1/ndiv)
#     Sheet["C"*string(ind)] = a
# end

