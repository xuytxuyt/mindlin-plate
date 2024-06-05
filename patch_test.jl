using ApproxOperator, JLD, XLSX

import BenchmarkExample: BenchmarkExample

include("import_patch_test.jl")
ndiv = 8
# elements, nodes = import_patch_test_fem("msh/patchtest_"*string(ndiv)*".msh");
elements, nodes = import_patch_test_fem("msh/patchtest_quad_"*string(ndiv)*".msh");

nₚ = length(nodes)

E = 1;
ν = 0.3;
h = 1
Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))

w(x,y) = -Dᵇ/Dˢ*8*x-Dᵇ/Dˢ*8*y+x^3+y^3+x^2*y+x*y^2
w₁(x,y) = -Dᵇ/Dˢ*8+3*x^2+2*x*y+y^2
w₂(x,y) = -Dᵇ/Dˢ*8+3*y^2+x^2+2*x*y
w₁₁(x,y) = 6*x+2*y
w₂₂(x,y) = 2*x+6*y
θ₁(x,y) = 3*x^2+2*x*y+y^2
θ₂(x,y) = 3*y^2+x^2+2*x*y
θ₁₁(x,y) = 6*x+2*y
θ₁₂(x,y) = 2*x+2*y
θ₂₂(x,y) = 2*x+6*y

# w(x,y) = x+y+x^2/2+x*y+y^2/2
# w₁(x,y) = 1+x+y
# w₂(x,y) = 1+x+y
# w₁₁(x,y) = 1
# w₂₂(x,y) = 1
# θ₁(x,y) = 1+x+y
# θ₂(x,y) = 1+x+y
# θ₁₁(x,y)  = 1
# θ₁₂(x,y)  = 1
# θ₂₂(x,y)  = 1

# n = 2
# w(x,y) = (x+y)^n
# w₁(x,y) = n*(x+y)^abs(n-1)
# w₂(x,y) = n*(x+y)^abs(n-1)
# w₁₁(x,y) = n*(n-1)*(x+y)^abs(n-2)
# w₂₂(x,y) = n*(n-1)*(x+y)^abs(n-2)
# m = 1
# θ₁(x,y) = (x+y)^m
# θ₂(x,y) = (x+y)^m
# θ₁₁(x,y)  = m*(x+y)^abs(m-1)
# θ₁₂(x,y)  = m*(x+y)^abs(m-1)
# θ₂₂(x,y)  = m*(x+y)^abs(m-1)
# θ₁₁₁(x,y) = m*(m-1)*(x+y)^abs(m-2)
# θ₁₁₂(x,y) = m*(m-1)*(x+y)^abs(m-2)
# θ₁₂₂(x,y) = m*(m-1)*(x+y)^abs(m-2)
# θ₁₂₁(x,y) = m*(m-1)*(x+y)^abs(m-2)
# θ₂₂₂(x,y) = m*(m-1)*(x+y)^abs(m-2)
# θ₂₂₁(x,y) = m*(m-1)*(x+y)^abs(m-2)

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
# ops[5](elements["Ω"],f)
ops[6](elements["Γ₁"],k,f)
ops[6](elements["Γ₂"],k,f)
ops[6](elements["Γ₃"],k,f)
ops[6](elements["Γ₄"],k,f)
# ops[7](elements["Γ₁"],k,f)
# ops[7](elements["Γ₂"],k,f)
# ops[7](elements["Γ₃"],k,f)
# ops[7](elements["Γ₄"],k,f)
# ops[8](elements["Γ₁"],f)
# ops[8](elements["Γ₂"],f)
# ops[8](elements["Γ₃"],f)
# ops[8](elements["Γ₄"],f)
ops[9](elements["Γ₁"],f)
ops[9](elements["Γ₂"],f)
ops[9](elements["Γ₃"],f)
ops[9](elements["Γ₄"],f)

d = (kᵇ+kˢ+k)\f
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d=>d₁)
# push!(nodes,:d=>d₂)
# push!(nodes,:d=>d₃)

set𝝭!(elements["Ωᵍ"])
set∇𝝭!(elements["Ωᵍ"])
prescribe!(elements["Ωᵍ"],:u=>(x,y,z)->w(x,y))
L₂ = ops[10](elements["Ωᵍ"])
a = log10(L₂)
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/patch_test.xlsx", mode="rw") do xf
#     Sheet = xf[5]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["F"*string(ind)] = log10(1/ndiv)
#     Sheet["G"*string(ind)] = L₂
#     Sheet["H"*string(ind)] = log10(L₂)
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