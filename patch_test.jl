using ApproxOperator, JLD,XLSX

import BenchmarkExample: BenchmarkExample

include("import_patch_test.jl")
ndiv = 3
elements, nodes = import_patch_test_fem("msh/patchtest_"*string(ndiv)*".msh");

nₚ = length(nodes)

E = 3E6;
ν = 0.3;
h = 1
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

n =2
w(x,y) = (x+y)^n
w₁(x,y) = n*(x+y)^abs(n-1)
w₂(x,y) = n*(x+y)^abs(n-1)
w₁₁(x,y) = n*(n-1)*(x+y)^abs(n-2)
w₂₂(x,y) = n*(n-1)*(x+y)^abs(n-2)
θ₁(x,y) = 1.0
θ₂(x,y) = 1.0
θ₁₁(x,y)  = 0.0
θ₁₂(x,y)  = 0.0
θ₂₂(x,y)  = 0.0
θ₁₁₁(x,y) = 0.0
θ₁₁₂(x,y) = 0.0
θ₁₂₂(x,y) = 0.0
θ₁₂₁(x,y) = 0.0
θ₂₂₂(x,y) = 0.0
θ₂₂₁(x,y) = 0.0
# θ₁(x,y) = n*(x+y)^abs(n-1)
# θ₂(x,y) = n*(x+y)^abs(n-1)
# θ₁₁(x,y) = n*(n-1)*(x+y)^abs(n-2)
# θ₁₂(x,y) = n*(n-1)*(x+y)^abs(n-2)
# θ₂₂(x,y) = n*(n-1)*(x+y)^abs(n-2)
# θ₁₁₁(x,y) = n*(n-1)*(n-2)*(x+y)^abs(n-3)
# θ₁₁₂(x,y) = n*(n-1)*(n-2)*(x+y)^abs(n-3)
# θ₁₂₂(x,y) = n*(n-1)*(n-2)*(x+y)^abs(n-3)
# θ₁₂₁(x,y) = n*(n-1)*(n-2)*(x+y)^abs(n-3)
# θ₂₂₂(x,y) = n*(n-1)*(n-2)*(x+y)^abs(n-3)
# θ₂₂₁(x,y) = n*(n-1)*(n-2)*(x+y)^abs(n-3)

# w(x,y) = (x^2+y^2)^2/64/Dᵇ-(x^2+y^2)*(1/Dˢ/4+1/32/Dᵇ)
# w₁(x,y) = x*(x^2+y^2)/16/Dᵇ-2*x*(1/Dˢ/4+1/32/Dᵇ)
# w₂(x,y) = y*(x^2+y^2)/16/Dᵇ-2*y*(1/Dˢ/4+1/32/Dᵇ)
# w₁₁(x,y) = (3*x^2+y^2)/16/Dᵇ-2*(1/Dˢ/4+1/32/Dᵇ)
# w₂₂(x,y) = (x^2+3*y^2)/16/Dᵇ-2*(1/Dˢ/4+1/32/Dᵇ)
# θ₁(x,y) = x*(x^2+y^2-1)/16/Dᵇ
# θ₂(x,y) = y*(x^2+y^2-1)/16/Dᵇ
# θ₁₁(x,y) = (3*x^2+y^2)/16/Dᵇ
# θ₁₂(x,y) = 2*x*y/16/Dᵇ
# θ₂₂(x,y) = (x^2+3*y^2)/16/Dᵇ

M₁₁(x,y)= -Dᵇ*(θ₁₁(x,y)+ν*θ₂₂(x,y))
M₁₂(x,y)= -Dᵇ*(1-ν)*θ₁₂(x,y)
M₂₂(x,y)= -Dᵇ*(ν*θ₁₁(x,y)+θ₂₂(x,y))
M₁₁₁(x,y)= -Dᵇ*(θ₁₁₁(x,y)+ν*θ₂₂₁(x,y))
M₁₂₂(x,y)= -Dᵇ*(1-ν)*θ₁₂₂(x,y)
M₁₂₁(x,y)= -Dᵇ*(1-ν)*θ₁₂₁(x,y)
M₂₂₂(x,y)= -Dᵇ*(ν*θ₁₁₂(x,y)+θ₂₂₂(x,y))

Q₁(x,y) = Dˢ*(w₁(x,y)-θ₁(x,y))
Q₂(x,y) = Dˢ*(w₂(x,y)-θ₂(x,y))
Q₁₁(x,y) = Dˢ*(w₁₁(x,y)-θ₁₁(x,y))
Q₂₂(x,y) = Dˢ*(w₂₂(x,y)-θ₂₂(x,y))

eval(prescribeForFem)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Γ₁"])
set𝝭!(elements["Γ₂"])
set𝝭!(elements["Γ₃"])
set𝝭!(elements["Γ₄"])

ops = [
    Operator{:∫κMγQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫γQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wqdΩ}(),
    Operator{:∫∇MQdΩ}(),
    Operator{:∫vwdΓ}(:α=>1e13*E),
    Operator{:∫vθdΓ}(:α=>1e13*E),
    Operator{:∫wVdΓ}(),
    Operator{:∫θMdΓ}(),
    Operator{:L₂}(:E=>E,:ν=>ν),
    Operator{:H₁}(:E=>E,:ν=>ν),

]
k = zeros(3*nₚ,3*nₚ)
kᵇ = zeros(3*nₚ,3*nₚ)
kˢ = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)
ops[2](elements["Ω"],kᵇ)
ops[3](elements["Ω"],kˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Ω"],f)
ops[6](elements["Γ₁"],k,f)
ops[6](elements["Γ₂"],k,f)
ops[6](elements["Γ₃"],k,f)
ops[6](elements["Γ₄"],k,f)
ops[7](elements["Γ₁"],k,f)
ops[7](elements["Γ₂"],k,f)
ops[7](elements["Γ₃"],k,f)
ops[7](elements["Γ₄"],k,f)
# ops[8](elements["Γ₁"],f)
# ops[8](elements["Γ₂"],f)
# ops[8](elements["Γ₃"],f)
# ops[8](elements["Γ₄"],f)
# ops[9](elements["Γ₁"],f)
# ops[9](elements["Γ₂"],f)
# ops[9](elements["Γ₃"],f)
# ops[9](elements["Γ₄"],f)

d = (kᵇ+kˢ+k)\f
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d=>d₁)
L₂ = ops[10](elements["Ω"])
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/patch_test.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndiv,index)+11
#     Sheet["B"*string(ind)] = log10(ndiv)
#     Sheet["C"*string(ind)] = L₂
#     Sheet["D"*string(ind)] = log10(L₂)
# end