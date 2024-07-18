using ApproxOperator, JLD, XLSX, LinearAlgebra

import BenchmarkExample: BenchmarkExample

include("import_SquarePlate.jl")
ndiv = 1
elements, nodes = import_SquarePlate_p("msh/SquarePlate_quad_"*string(ndiv)*".msh");
nₚ = 28
E = 1
ν = 1
h = 1
L = 1

Dᵇ = E*h^3/12/(1-ν^2)
push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])

ops = [
    Operator{:∫κMγQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫γQdΩ}(:E=>E,:ν=>ν,:h=>h),
]
k = zeros(3*nₚ,3*nₚ)
kᵇ = zeros(3*nₚ,3*nₚ)
kˢ = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)

ops[2](elements["Ω"],kᵇ) 
ops[3](elements["Ω"],kˢ)
rank(kˢ)
# for i in (1:3*nₚ)
#     for j in (1:3*nₚ)
#         XLSX.openxlsx("./xlsx/rank.xlsx", mode="rw") do xf
#            Sheet = xf[6]
#            Sheet[i,j] = kˢ[i,j]
#         end
#     end
# end

