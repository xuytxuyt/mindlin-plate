using ApproxOperator, JLD, XLSX, Printf
using CairoMakie
# using SparseArrays, Pardiso
import BenchmarkExample: BenchmarkExample
include("import_SquarePlate.jl")
ndiv   = 8
indices = 5:5

n_eig_nonzeros = zeros(Int,length(indices))
n_eig_real = zeros(Int,length(indices))
min_eig_nonzeros = zeros(length(indices))
min_eig_real = zeros(length(indices))

for (i,n) in enumerate(indices)

elements, nodes, nodes_s = import_infsup_linear_mix("msh/SquarePlate/SquarePlate_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_"*string(n)*".msh",n);

# elements, nodes, nodes_s, Ω, sp, type = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_quad_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_quad_"*string(ndivs)*".msh");

# elements, nodes, nodes_s, Ω, sp, type = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_tri6_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_tri6_"*string(ndivs)*".msh");


# elements, nodes, nodes_s, Ω, sp, type = import_SquarePlate_mix("msh/SquarePlate/SquarePlate_quad8_"*string(ndiv)*".msh","msh/SquarePlate/SquarePlate_quad8_"*string(ndivs)*".msh");


nᵇ = length(nodes)
nˢ = length(nodes_s)
nₑ = length(elements["Ω"])
E = 10.92e6
ν = 0.3
h = 0.001
L = 1.0
# ps = MKLPardisoSolver()

Dᵇ = E*h^3/12/(1-ν^2)
Dˢ = 5/6*E*h/(2*(1+ν))
w(x,y) = 1/3*x^3*(x-1)^3*y^3*(y-1)^3-2*h^2/(5*(1-ν))*(y^3*(y-1)^3*x*(x-1)*(5*x^2-5*x+1)+x^3*(x-1)^3*y*(y-1)*(5*y^2-5*y+1))
θ₁(x,y) = y^3*(y-1)^3*x^2*(x-1)^2*(2*x-1)
θ₂(x,y) = x^3*(x-1)^3*y^2*(y-1)^2*(2*y-1)
F(x,y) = E*h^3/(12*(1-ν^2))*(12*y*(y-1)*(5*x^2-5*x+1)*(2*y^2*(y-1)^2+x*(x-1)*(5*y^2-5*y+1))+12*x*(x-1)*(5*y^2-5*y+1)*(2*x^2*(x-1)^2+y*(y-1)*(5*x^2-5*x+1)))

w₁(x,y) = (x-1)^2*x^2*(2*x-1)*(y-1)^3*y^3-2*h^2/(5*(1-ν))*((20*x^3-30*x^2+12*x-1)*(y-1)^3*y^3+3*(x-1)^2*x^2*(2*x-1)*(y-1)*y*(5*y^2-5*y+1))
w₂(x,y) = (x-1)^3*x^3*(y-1)^2*y^2*(2*y-1)-2*h^2/(5*(1-ν))*(3*(x-1)*x*(5*x^2-5*x+1)*(y-1)^2*y^2*(2*y-1)+x^3*(x-1)^3*(20*y^3-30*y^2+12*y-1))
θ₁₁(x,y) = 2*(x-1)*x*(5*x^2-5*x+1)*(y-1)^3*y^3
θ₁₂(x,y) = 3*(x-1)^2*x^2*(2*x-1)*(y-1)^2*y^2*(2*y-1)
θ₂₂(x,y) = 2*(x-1)^3*x^3*(y-1)*y*(5*y^2-5*y+1)
M₁₁(x,y)= -Dᵇ*(θ₁₁(x,y)+ν*θ₂₂(x,y))
M₁₂(x,y)= -Dᵇ*(1-ν)*θ₁₂(x,y)
M₂₂(x,y)= -Dᵇ*(ν*θ₁₁(x,y)+θ₂₂(x,y))
Q₁(x,y) = Dˢ*(w₁(x,y)-θ₁(x,y))
Q₂(x,y) = Dˢ*(w₂(x,y)-θ₂(x,y))


set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γˡ"])
set𝝭!(elements["Γʳ"])

prescribe!(elements["Γᵇ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γᵗ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γˡ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γʳ"],:g=>(x,y,z)->w(x,y))
prescribe!(elements["Γᵇ"],:θ₁=>(x,y,z)->θ₁(x,y))
prescribe!(elements["Γᵗ"],:θ₁=>(x,y,z)->θ₁(x,y))
prescribe!(elements["Γˡ"],:θ₁=>(x,y,z)->θ₁(x,y))
prescribe!(elements["Γʳ"],:θ₁=>(x,y,z)->θ₁(x,y))
prescribe!(elements["Γᵇ"],:θ₂=>(x,y,z)->θ₂(x,y))
prescribe!(elements["Γᵗ"],:θ₂=>(x,y,z)->θ₂(x,y))
prescribe!(elements["Γˡ"],:θ₂=>(x,y,z)->θ₂(x,y))
prescribe!(elements["Γʳ"],:θ₂=>(x,y,z)->θ₂(x,y))
prescribe!(elements["Ω"],:q=>(x,y,z)->F(x,y))

ops = [
    Operator{:∫κMdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wQdΩ}(),
    Operator{:∫QQdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫wqdΩ}(),
    Operator{:∫vwdΓ}(:α=>1e13*E),
    Operator{:∫vθ₁dΓ}(:α=>1e13*E),
    Operator{:∫vθ₂dΓ}(:α=>1e13*E),
    Operator{:L₂_ThickPlate}(:E=>E,:ν=>ν),
    Operator{:L₂_ThickPlate_Q}(:E=>E,:ν=>ν),
    Operator{:∫θM₁dΓ}(),
    Operator{:∫θM₂dΓ}(),
    Operator{:∫wVdΓ}(),
]
kᵇ = zeros(3*nᵇ,3*nᵇ)
kʷˢ = zeros(3*nᵇ,2*nˢ)
kˢˢ = zeros(2*nˢ,2*nˢ)
f = zeros(3*nᵇ)
# d = zeros(3*nᵇ+2*nˢ)

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


# k = kʷˢ*inv(kˢˢ)*kʷˢ'
val = eigvals(kʷˢ*(kˢˢ\kʷˢ'))
println(val)
# println(log10(a[3*nᵇ-2nˢ+1]))
# println(a[3*nᵇ-2nˢ+1])

val_sign = zeros(3*nᵇ)
for (ii,v) in enumerate(val)
    if v isa Real
        val_sign[ii] = sign(v)
    else
        val_sign[ii] = sign(v.re) < -1e-8 ? -1.0 : 1.0
    end
end
val_real = val_sign .* abs.(val)
val_abs = abs.(val)
# println("Sorted Eigenvalue")
val_sort = sort(val_abs)
# println.(val_sort[2*nᵤ-nₚ.+(-2:4)]);

n_eig_real[i] = count(x-> abs(x)>1e-8, val_real)
n_eig_nonzeros[i] = count(x-> x > 1e-8,val_sort)
min_eig_real[i] = min(val_real[val_real.>1e-8]...)
min_eig_nonzeros[i] = val_sort[3*nᵇ - n_eig_nonzeros[i] + 1]

end

XLSX.openxlsx("./xlsx/infsup.xlsx", mode = "rw") do xf
    sheet = xf[1]
    for (n,n_eig_r,min_eig_r,n_eig_n,min_eig_n) in zip(indices,n_eig_real,min_eig_real,n_eig_nonzeros,min_eig_nonzeros)
        sheet["A"*string(n)] = n
        sheet["B"*string(n)] = n_eig_r
        sheet["C"*string(n)] = min_eig_r^0.5
        sheet["D"*string(n)] = n_eig_n
        sheet["E"*string(n)] = min_eig_n^0.5
    end
end