using ApproxOperator, JLD, XLSX

import BenchmarkExample: BenchmarkExample

include("import_Circular.jl")

ndiv  = 8
ndivs = 8
elements, nodes, nodes_s= import_Circular_mix("msh/circular_"*string(ndiv)*".msh","msh/circular_"*string(ndivs)*".msh");
nᵇ = length(nodes)
nˢ = length(nodes_s)
E = BenchmarkExample.Circular.𝐸
ν = BenchmarkExample.Circular.𝜈
h = BenchmarkExample.Circular.ℎ
F = BenchmarkExample.Circular.𝐹
w = BenchmarkExample.Circular.𝑤
R = BenchmarkExample.Circular.𝑅
θ₁ = BenchmarkExample.Circular.𝜃₁
θ₂ = BenchmarkExample.Circular.𝜃₂
𝑣 = BenchmarkExample.Circular.𝑣ᴱ
eval(prescribeForSSUniformLoading)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γᵉ"])
set𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])

ops = [
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wQdΩ}(),
    Operator{:∫QQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wqdΩ}(),
    Operator{:∫vwdΓ}(:α=>1e13*E),
    Operator{:∫vθ₁dΓ}(:α=>1e13*E),
    Operator{:∫vθ₂dΓ}(:α=>1e13*E),
]
kᵇ = zeros(3*nᵇ,3*nᵇ)
kʷˢ = zeros(3*nᵇ,2*nˢ)
kˢˢ = zeros(2*nˢ,2*nˢ)
f = zeros(3*nᵇ)
ops[1](elements["Ω"],kᵇ)
ops[2](elements["Ω"],elements["Ωˢ"],kʷˢ)
ops[3](elements["Ωˢ"],kˢˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Γᵉ"],kᵇ,f)
ops[6](elements["Γˡ"],kᵇ,f)
ops[7](elements["Γᵇ"],kᵇ,f)

ops𝐴 = Operator{:SphericalShell_𝐴}()
k = [kᵇ kʷˢ;kʷˢ' kˢˢ]
f = [f;zeros(2*nˢ)]

d = k\f
d₁ = d[1:3:3*nᵇ]
d₂ = d[2:3:3*nᵇ]
d₃ = d[3:3:3*nᵇ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = ops𝐴(elements["𝐴"])
# e = abs(w[1]-𝑣)
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/Circular_SSUL.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndivs,index)+1
#     Sheet["F"*string(ind)] = log10(5/ndiv)
#     Sheet["G"*string(ind)] = w
#     Sheet["H"*string(ind)] = log10(e)
# end