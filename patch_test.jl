using ApproxOperator, JLD

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

n =3
w(x,y) = -Dᵇ/Dˢ*12*x-Dᵇ/Dˢ*12*y+(x+y)^n
w₁(x,y) = -Dᵇ/Dˢ*12+n*(x+y)^abs(n-1)
w₂(x,y) = -Dᵇ/Dˢ*12+n*(x+y)^abs(n-1)
w₁₁(x,y) = n*(n-1)*(x+y)^abs(n-2)
w₂₂(x,y) = n*(n-1)*(x+y)^abs(n-2)
θ₁(x,y) = n*(x+y)^abs(n-1)
θ₂(x,y) = n*(x+y)^abs(n-1)
θ₁₁(x,y) = n*(n-1)*(x+y)^abs(n-2)
θ₁₂(x,y) = n*(n-1)*(x+y)^abs(n-2)
θ₂₂(x,y) = n*(n-1)*(x+y)^abs(n-2)

M₁₁(x,y)= -Dᵇ*(θ₁₁(x,y)+ν*θ₂₂(x,y))
M₁₂(x,y)= -Dᵇ*(1-ν)*θ₁₂(x,y)
M₂₂(x,y)= -Dᵇ*(ν*θ₁₁(x,y)+θ₂₂(x,y))

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
    Operator{:∫vwdΓ}(:α=>1e8*E),
    Operator{:∫vθdΓ}(:α=>1e8*E),
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
ops[5](elements["Γ₁"],k,f)
ops[5](elements["Γ₂"],k,f)
ops[5](elements["Γ₃"],k,f)
ops[5](elements["Γ₄"],k,f)
ops[6](elements["Γ₁"],k,f)
ops[6](elements["Γ₂"],k,f)
ops[6](elements["Γ₃"],k,f)
ops[6](elements["Γ₄"],k,f)
ops[7](elements["Γ₁"],f)
ops[7](elements["Γ₂"],f)
ops[7](elements["Γ₃"],f)
ops[7](elements["Γ₄"],f)
ops[8](elements["Γ₁"],f)
ops[8](elements["Γ₂"],f)
ops[8](elements["Γ₃"],f)
ops[8](elements["Γ₄"],f)

d = (kᵇ+kˢ+k)\f
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d=>d₂)
L₂ = ops[9](elements["Ω"])