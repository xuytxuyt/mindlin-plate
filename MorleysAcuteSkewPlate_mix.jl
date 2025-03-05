using ApproxOperator, JLD, XLSX, Printf
import BenchmarkExample: BenchmarkExample
using SparseArrays, Pardiso
using CairoMakie

include("import_MorleysAcuteSkewPlate.jl")
ndiv  = 16
ndivs = 12
# elements, nodes, nodes_s, Ω, sp, type= import_MorleysAcuteSkewPlate_mix("msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_"*string(ndiv)*".msh","msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω, sp, type= import_MorleysAcuteSkewPlate_mix("msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad_"*string(ndiv)*".msh","msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω, sp, type= import_MorleysAcuteSkewPlate_mix("msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_tri6_"*string(ndiv)*".msh","msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_tri6_"*string(ndivs)*".msh");
elements, nodes, nodes_s, Ω, sp, type= import_MorleysAcuteSkewPlate_mix("msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad8_"*string(ndiv)*".msh","msh/MorleysAcuteSkewPlate/MorleysAcuteSkewPlate_quad8_"*string(ndivs)*".msh");

nᵇ = length(nodes)
nˢ = length(nodes_s)

E = 1085.0
ν = 0.31
h = 0.1
F = 1.0
w = 0.0
L = 100.0
θ₁ = 0.0
θ₂ = 0.0
𝑣 = 0.7945
Dᵇ = E*h^3/12/(1-ν^2)
# ps = MKLPardisoSolver()

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
kᵇ = spzeros(3*nᵇ,3*nᵇ)
kʷˢ = spzeros(3*nᵇ,2*nˢ)
kˢˢ = spzeros(2*nˢ,2*nˢ)
f = zeros(3*nᵇ)

ops[1](elements["Ω"],kᵇ)
ops[2](elements["Ω"],elements["Ωˢ"],kʷˢ)
ops[3](elements["Ωˢ"],kˢˢ)
ops[4](elements["Ω"],f)
ops[5](elements["Γᵇ"],kᵇ,f)
ops[5](elements["Γᵗ"],kᵇ,f)
ops[6](elements["Γᵇ"],kᵇ,f)
ops[6](elements["Γᵗ"],kᵇ,f)


ops𝐴 = Operator{:SphericalShell_𝐴}()
k = [kᵇ kʷˢ;kʷˢ' kˢˢ]
f = [f;zeros(2*nˢ)]

d = k\f
d₁ = d[1:3:3*nᵇ]
d₂ = d[2:3:3*nᵇ]
d₃ = d[3:3:3*nᵇ]
s₁ = d[3*nᵇ+1:2:3*nᵇ+2*nˢ]
s₂ = d[3*nᵇ+2:2:3*nᵇ+2*nˢ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
push!(nodes_s,:q₁=>s₁,:q₂=>s₂)

w = ops𝐴(elements["𝐴"])
wᶜ= w*10^2*Dᵇ/(F*L^4)

# println(wᶜ)
# index = 20:64
# XLSX.openxlsx("./xlsx/SquarePlate.xlsx", mode="rw") do xf
#     Sheet = xf[2]
#     ind = findfirst(n->n==ndivs,index)+1
#     Sheet["A"*string(ind)] = nˢ
#     Sheet["B"*string(ind)] = log10(abs(1-abs(wᶜ[1]/𝑣)))
# end

# println(wᶜ)
# index = [2,4,6,8,16,24,32,48,64]
# XLSX.openxlsx("./xlsx/SquarePlate.xlsx", mode="rw") do xf
#     Sheet = xf[3]
#     ind = findfirst(n->n==ndiv,index)+1
#     Sheet["A"*string(ind)] = ndiv
#     Sheet["B"*string(ind)] = abs(wᶜ[1]/𝑣)
# end



ind = 100

# xs = LinRange(0, 1, ind)
# ys = LinRange(0, 1, ind)
xs₁ = zeros(ind,ind)
ys₁ = zeros(ind,ind)
zs₁ = zeros(ind,ind)
xs₂ = zeros(ind,ind)
ys₂ = zeros(ind,ind)
zs₂ = zeros(ind,ind)
xs₃ = zeros(ind,ind)
ys₃ = zeros(ind,ind)
zs₃ = zeros(ind,ind)
𝗠 = zeros(21)

for i in 1:ind
    for j in 1:ind
        Δy = 50*3^0.5/(ind-1)
        Δx = (1/3^0.5*(i-1)Δy)/(ind-1)
        x₁ = 50+(j-(ind))*Δx
        y₁ = 50*3^0.5-(i-1)*Δy
        indices = sp(x₁,y₁,0.0)
        ni = length(indices)
        𝓒 = [nodes_s[i] for i in indices]
        data = Dict([:x=>(2,[x₁]),:y=>(2,[y₁]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
        ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
        𝓖 = [ξ]
        a = type(𝓒,𝓖)
        set𝝭!(a)
        q = 0.0
        N = ξ[:𝝭]
        for (k,xₖ) in enumerate(𝓒)
            q += N[k]*xₖ.q₁
            # q += N[k]*xₖ.q₂
        end
        xs₁[i,j] = x₁
        ys₁[i,j] = y₁
        zs₁[i,j] = q
    end
 end

for (i,x₁) in enumerate(LinRange(50,100, ind))
    for (j,y₁) in enumerate(LinRange(0.0, 50*3^0.5, ind))
        indices = sp(x₁,y₁,0.0)
        ni = length(indices)
        𝓒 = [nodes_s[i] for i in indices]
        data = Dict([:x=>(2,[x₁]),:y=>(2,[y₁]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
        ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
        𝓖 = [ξ]
        a = type(𝓒,𝓖)
        set𝝭!(a)
        q = 0.0
        N = ξ[:𝝭]
        for (k,xₖ) in enumerate(𝓒)
            q += N[k]*xₖ.q₁
            # q += N[k]*xₖ.q₂
        end
        xs₂[i,j] = x₁
        ys₂[i,j] = y₁
        zs₂[i,j] = q
    end
 end


 for i in 1:ind
    for j in 1:ind
        Δy = 50*3^0.5/(ind-1)
        Δx = (1/3^0.5*(ind-i)Δy)/(ind-1)
        x₁ = 100+((ind)-j)*Δx
        y₁ = 50*3^0.5-(i-1)*Δy
        indices = sp(x₁,y₁,0.0)
        ni = length(indices)
        𝓒 = [nodes_s[i] for i in indices]
        data = Dict([:x=>(2,[x₁]),:y=>(2,[y₁]),:z=>(2,[0.0]),:𝝭=>(4,zeros(ni)),:𝗠=>(0,𝗠)])
        ξ = 𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0), data)
        𝓖 = [ξ]
        a = type(𝓒,𝓖)
        set𝝭!(a)
        q = 0.0
        N = ξ[:𝝭]
        for (k,xₖ) in enumerate(𝓒)
            q += N[k]*xₖ.q₁
            # q += N[k]*xₖ.q₂
        end
        xs₃[i,j] = x₁
        ys₃[i,j] = y₁
        zs₃[i,j] = q
    end
 end

fig = Figure()
ax = Axis(fig[1,1], 
    aspect = DataAspect(), 
    xticksvisible = false,
    xticklabelsvisible=false, 
    yticksvisible = false, 
    yticklabelsvisible=false,
)
hidespines!(ax)
hidedecorations!(ax)
# 

surface!(xs₁,ys₁,zeros(ind,ind),color=zs₁,colorrange=(-200,200),shading = false,colormap=:rainbow)
# contour!(zs₁,levels=-1360:300:1360,color=:azure)
surface!(xs₂,ys₂,zeros(ind,ind),color=zs₂,colorrange=(-200,200),shading = false,colormap=:rainbow)
# contour!(zs₂,levels=-1360:300:1360,color=:azure)
surface!(xs₃,ys₃,zeros(ind,ind),color=zs₃,colorrange=(-200,200),shading = false,colormap=:rainbow)
# contour!(zs₃,levels=-1360:200:1360,color=:azure)
# Colorbar(fig[1,2], limits=(-1360,1360), colormap=:lightrainbow)
# save("./png/MorleysAcuteSkewPlate_tri3_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/MorleysAcuteSkewPlate_colorbar.png",fig, px_per_unit = 3.0)
# save("./png/MorleysAcuteSkewPlate_tri3_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/MorleysAcuteSkewPlate_tri6_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/MorleysAcuteSkewPlate_tri6_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/MorleysAcuteSkewPlate_quad4_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/MorleysAcuteSkewPlate_quad4_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)
# save("./png/MorleysAcuteSkewPlate_quad8_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/MorleysAcuteSkewPlate_quad8_q2_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 10.0)

fig