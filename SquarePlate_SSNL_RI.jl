using ApproxOperator, JLD, XLSX

import BenchmarkExample: BenchmarkExample

include("import_SquarePlate.jl")
ndiv = 60
elements, nodes, nodes_s = import_SquarePlate_quad_RI("msh/SquarePlate_quad_"*string(ndiv)*".msh","msh/SquarePlate_quad_"*string(ndiv)*".msh");
nᵇ = length(nodes)
nˢ = length(nodes_s)
E = BenchmarkExample.SquarePlate.𝐸
ν = BenchmarkExample.SquarePlate.𝜈
h = BenchmarkExample.SquarePlate.ℎ
L = BenchmarkExample.SquarePlate.𝐿

Dᵇ = E*h^3/12/(1-ν^2)
w(x,y) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
θ₁(x,y) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
θ₂(x,y) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
F(x,y) = E*h^3/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))
eval(prescribeForSSNonUniformLoading)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γˡ"])
set𝝭!(elements["Γʳ"])

ops = [
    Operator{:∫κMγQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫γQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wqdΩ}(),
    Operator{:∫vwdΓ}(:α=>1e13*E),
    Operator{:∫vθ₁dΓ}(:α=>1e13*E),
    Operator{:∫vθ₂dΓ}(:α=>1e13*E),
    Operator{:L₂}(:E=>E,:ν=>ν),
]
k = zeros(3*nᵇ,3*nᵇ)
kᵇ = zeros(3*nᵇ,3*nᵇ)
kˢ = zeros(3*nˢ,3*nˢ)
f = zeros(3*nᵇ)
# ops[1](elements["Ω"],k)
ops[2](elements["Ω"],kᵇ)
ops[3](elements["Ωˢ"],kˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Γᵇ"],k,f)
ops[5](elements["Γᵗ"],k,f)
ops[5](elements["Γˡ"],k,f)
ops[5](elements["Γʳ"],k,f)
ops[6](elements["Γᵇ"],k,f)
ops[6](elements["Γᵗ"],k,f)
ops[6](elements["Γˡ"],k,f)
ops[6](elements["Γʳ"],k,f)
ops[7](elements["Γᵇ"],k,f)
ops[7](elements["Γᵗ"],k,f)
ops[7](elements["Γˡ"],k,f)
ops[7](elements["Γʳ"],k,f)

d = (kᵇ+kˢ+k)\f
d₁ = d[1:3:3*nᵇ]
d₂ = d[2:3:3*nᵇ]
d₃ = d[3:3:3*nᵇ]

push!(nodes,:d=>d₁)
set𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωᵍ"])
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->w(x,y))
L₂ = ops[8](elements["Ωᵍ"])
a = log10(L₂)
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/SquarePlate_UniformLoading.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndiv,index)+11
#     Sheet["B"*string(ind)] = log10(1/ndiv)
#     Sheet["C"*string(ind)] = a
# end