using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample

include("import_patch_test.jl")
ndiv = 8
elements, nodes = import_patch_test_fem("msh/patchtest_"*string(ndiv)*".msh");

nₚ = length(nodes)

E = 1;
ν = 0.3;
h = 0.1
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

n = 4
w(x,y) = (x+y)^n
w₁(x,y) = n*(x+y)^abs(n-1)
w₂(x,y) = n*(x+y)^abs(n-1)
w₁₁(x,y) = n*(n-1)*(x+y)^abs(n-2)
w₂₂(x,y) = n*(n-1)*(x+y)^abs(n-2)
m = 2
θ₁(x,y) = (x+y)^m
θ₂(x,y) = (x+y)^m
θ₁₁(x,y)  = m*(x+y)^abs(m-1)
θ₁₂(x,y)  = m*(x+y)^abs(m-1)
θ₂₂(x,y)  = m*(x+y)^abs(m-1)
θ₁₁₁(x,y) = m*(m-1)*(x+y)^abs(m-2)
θ₁₁₂(x,y) = m*(m-1)*(x+y)^abs(m-2)
θ₁₂₂(x,y) = m*(m-1)*(x+y)^abs(m-2)
θ₁₂₁(x,y) = m*(m-1)*(x+y)^abs(m-2)
θ₂₂₂(x,y) = m*(m-1)*(x+y)^abs(m-2)
θ₂₂₁(x,y) = m*(m-1)*(x+y)^abs(m-2)

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
# ops[1](elements["Ω"],k)
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
# push!(nodes,:d=>d₂)
# push!(nodes,:d=>d₃)
L₂ = ops[10](elements["Ω"])
a = log10(L₂)
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/patch_test.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["J"*string(ind)] = log10(1/ndiv)
#     Sheet["K"*string(ind)] = L₂
#     Sheet["L"*string(ind)] = log10(L₂)
# end

# d = zeros(3*nₚ)
# for (i,node) in enumerate(nodes)
#     x = node.x
#     y = node.y
#     z = node.z

#     d[3*i-2] = w(x,y)
#     d[3*i-1] = θ₁(x,y)
#     d[3*i] = θ₂(x,y)
# end

# k
# kᵇ*d - f
# kˢ*d - f
# err1 = kᵇ*d