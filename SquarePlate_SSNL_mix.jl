using ApproxOperator, JLD, XLSX

import BenchmarkExample: BenchmarkExample
include("import_SquarePlate.jl")
ndiv  = 20
ndivs = 20
elements, nodes, nodes_s= import_SquarePlate_mix("msh/SquarePlate_"*string(ndiv)*".msh","msh/SquarePlate_"*string(ndivs)*".msh");
# elements, nodes, nodes_s= import_SquarePlate_mix("msh/SquarePlate_bubble_"*string(ndiv)*".msh","msh/SquarePlate_bubble_"*string(ndivs)*".msh");
nᵇ = length(nodes)
nˢ = length(nodes_s)

E = BenchmarkExample.SquarePlate.𝐸
ν = BenchmarkExample.SquarePlate.𝜈
h = BenchmarkExample.SquarePlate.ℎ
L = BenchmarkExample.SquarePlate.𝐿
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
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wQdΩ}(),
    Operator{:∫QQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wqdΩ}(),
    Operator{:∫vwdΓ}(:α=>1e13*E),
    Operator{:∫vθ₁dΓ}(:α=>1e13*E),
    Operator{:∫vθ₂dΓ}(:α=>1e13*E),
    Operator{:L₂_ThickPlate}(:E=>E,:ν=>ν),
    Operator{:L₂}(:E=>E,:ν=>ν),
]
kᵇ = zeros(3*nᵇ,3*nᵇ)
kʷˢ = zeros(3*nᵇ,2*nˢ)
kˢˢ = zeros(2*nˢ,2*nˢ)
f = zeros(3*nᵇ)

ops[1](elements["Ω"],kᵇ)
ops[2](elements["Ω"],elements["Ωˢ"],kʷˢ)
ops[3](elements["Ωˢ"],kˢˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Γᵇ"],kᵇ,f)
ops[5](elements["Γᵗ"],kᵇ,f)
ops[5](elements["Γˡ"],kᵇ,f)
ops[5](elements["Γʳ"],kᵇ,f)
ops[6](elements["Γᵇ"],kᵇ,f)
ops[6](elements["Γᵗ"],kᵇ,f)
ops[6](elements["Γˡ"],kᵇ,f)
ops[6](elements["Γʳ"],kᵇ,f)
ops[7](elements["Γᵇ"],kᵇ,f)
ops[7](elements["Γᵗ"],kᵇ,f)
ops[7](elements["Γˡ"],kᵇ,f)
ops[7](elements["Γʳ"],kᵇ,f)

k = [kᵇ kʷˢ;kʷˢ' kˢˢ]
f = [f;zeros(2*nˢ)]

d = k\f
d₁ = d[1:3:3*nᵇ]
d₂ = d[2:3:3*nᵇ]
d₃ = d[3:3:3*nᵇ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
set𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωᵍ"])
prescribe!(elements["Ωᵍ"],:w=>(x,y,z)->w(x,y))
prescribe!(elements["Ωᵍ"],:θ₁=>(x,y,z)->θ₁(x,y))
prescribe!(elements["Ωᵍ"],:θ₂=>(x,y,z)->θ₂(x,y))
L₂ = ops[8](elements["Ωᵍ"])
a = log10(L₂)
# println(wᶜ)
# e = abs(wᶜ[1]-𝑣)
# index = [270,280,290,300,310,320,330,340,350,360,370,380,390,400]
# XLSX.openxlsx("./xlsx/SquarePlate_UniformLoading.xlsx", mode="rw") do xf
#     Sheet = xf[4]
#     ind = findfirst(n->n==ndivs,index)+1
#     Sheet["H"*string(ind)] = nˢ
#     Sheet["I"*string(ind)] = a
# end