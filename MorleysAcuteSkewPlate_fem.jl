using ApproxOperator, JLD, XLSX, LinearAlgebra

import BenchmarkExample: BenchmarkExample

include("import_MorleysAcuteSkewPlate.jl")
ndiv = 2
elements, nodes = import_MorleysAcuteSkewPlate("msh/MorleysAcuteSkewPlate_"*string(ndiv)*".msh");
nₚ = length(nodes)

E = BenchmarkExample.MorleysAcuteSkewPlate.𝐸
ν = BenchmarkExample.MorleysAcuteSkewPlate.𝜈
h = BenchmarkExample.MorleysAcuteSkewPlate.ℎ
F = BenchmarkExample.MorleysAcuteSkewPlate.𝐹
w = BenchmarkExample.MorleysAcuteSkewPlate.𝑤
L = BenchmarkExample.MorleysAcuteSkewPlate.𝐿
θ₁ = BenchmarkExample.MorleysAcuteSkewPlate.𝜃₁
θ₂ = BenchmarkExample.MorleysAcuteSkewPlate.𝜃₂
𝑣 = BenchmarkExample.MorleysAcuteSkewPlate.𝑣ᴱ
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
ops[5](elements["Γᵗ"],k,f)
ops[6](elements["Γᵇ"],k,f)
ops[6](elements["Γᵗ"],k,f)

ops𝐴 = Operator{:SphericalShell_𝐴}()

d = (kᵇ+kˢ+k)\f
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = ops𝐴(elements["𝐴"])
wᶜ= w*10^2*Dᵇ/(F*L^4)

# println(wᶜ)
index = [2,4,6,8,16,24,32,48,64]
XLSX.openxlsx("./xlsx/SquarePlate.xlsx", mode="rw") do xf
    Sheet = xf[3]
    ind = findfirst(n->n==ndiv,index)+1
    Sheet["A"*string(ind)] = ndiv
    Sheet["B"*string(ind)] = abs(wᶜ[1]/𝑣)
end