using ApproxOperator, JLD, XLSX

import BenchmarkExample: BenchmarkExample

include("import_SquarePlate.jl")
ndiv = 8
elements, nodes = import_SquarePlate("msh/QuarterSquarePlate_"*string(ndiv)*".msh");
nₚ = length(nodes)

E = BenchmarkExample.SquarePlate.𝐸
ν = BenchmarkExample.SquarePlate.𝜈
h = BenchmarkExample.SquarePlate.ℎ
L = BenchmarkExample.SquarePlate.𝐿
F = BenchmarkExample.SquarePlate.𝐹

Dᵇ = E*h^3/12/(1-ν^2)
eval(prescribeForSSUniformLoading)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γˡ"])
set𝝭!(elements["Γʳ"])
set𝝭!(elements["𝐴"])

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

k = zeros(3*nₚ,3*nₚ)
kᵇ = zeros(3*nₚ,3*nₚ)
kˢ = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)
# ops[1](elements["Ω"],k)
ops[2](elements["Ω"],kᵇ)
ops[3](elements["Ω"],kˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Γᵇ"],k,f)
ops[5](elements["Γˡ"],k,f)
ops[6](elements["Γᵇ"],k,f)
ops[6](elements["Γʳ"],k,f)
ops[7](elements["Γᵗ"],k,f)
ops[7](elements["Γˡ"],k,f)

d = (kᵇ+kˢ+k)\f
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁)
ops𝐴 = Operator{:SphericalShell_𝐴}()

w = ops𝐴(elements["𝐴"])
wᶜ= w*10^3*Dᵇ/(F*L^4)
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/SquarePlate_UniformLoading.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndiv,index)+6
#     Sheet["B"*string(ind)] = log10(1/ndiv)
#     Sheet["C"*string(ind)] = a
# end