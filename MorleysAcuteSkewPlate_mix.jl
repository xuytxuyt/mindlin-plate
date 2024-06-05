using ApproxOperator, JLD, XLSX

import BenchmarkExample: BenchmarkExample

include("import_MorleysAcuteSkewPlate.jl")
ndiv = 8
elements, nodes, nodes_s= import_MorleysAcuteSkewPlate_mix("msh/MorleysAcuteSkewPlate_"*string(ndiv)*".msh","msh/MorleysAcuteSkewPlate_"*string(ndiv)*".msh");
nᵇ = length(nodes)
nˢ = length(nodes_s)

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
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γˡ"])
set𝝭!(elements["Γʳ"])
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
ops[5](elements["Γᵇ"],kᵇ,f)
ops[5](elements["Γᵗ"],kᵇ,f)
ops[5](elements["Γˡ"],kᵇ,f)
ops[5](elements["Γʳ"],kᵇ,f)


ops𝐴 = Operator{:SphericalShell_𝐴}()
k = [kᵇ kʷˢ;kʷˢ' kˢˢ]
f = [f;zeros(2*nˢ)]

d = k\f
d₁ = d[1:3:3*nᵇ]
d₂ = d[2:3:3*nᵇ]
d₃ = d[3:3:3*nᵇ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = ops𝐴(elements["𝐴"])
wᶜ= w*10^3*Dᵇ/(F*L^4)
# println(wᶜ)
e = abs(wᶜ[1]-𝑣)
index = [8,16,32,64]
XLSX.openxlsx("./xlsx/MorleysAcuteSkewPlate.xlsx", mode="rw") do xf
    Sheet = xf[1]
    ind = findfirst(n->n==ndiv,index)+11
    Sheet["F"*string(ind)] = log10(100/ndiv)
    Sheet["G"*string(ind)] = wᶜ
    Sheet["H"*string(ind)] = log10(e)
end