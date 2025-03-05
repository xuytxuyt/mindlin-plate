using ApproxOperator, JLD, XLSX, GLMakie 

import BenchmarkExample: BenchmarkExample
using SparseArrays, Pardiso

include("import_Circular.jl")
ndiv  = 16
ndivs = 16
elements, nodes, nodes_s, Ω, sp, type= import_Circular_mix("msh/circular/circular_"*string(ndiv)*".msh","msh/circular/circular_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω, sp, type= import_Circular_mix("msh/circular/circular_quad_"*string(ndiv)*".msh","msh/circular/circular_quad_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω, sp, type= import_Circular_mix("msh/circular/circular_tri6_"*string(ndiv)*".msh","msh/circular/circular_tri6_"*string(ndivs)*".msh");
# elements, nodes, nodes_s, Ω, sp, type= import_Circular_mix("msh/circular/circular_quad8_"*string(ndiv)*".msh","msh/circular/circular_quad8_"*string(ndivs)*".msh");
nᵇ = length(nodes)
nˢ = length(nodes_s)
E = 10.92
ν = 0.3
h = 0.1
F = 1.0
w = 0.0
R = 5.0
θ₁ = 0.0
θ₂ = 0.0
𝑣 = 39831.0
eval(prescribeForSSUniformLoading)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set𝝭!(elements["Ωˢ"])
set∇𝝭!(elements["Ωˢ"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γᵉ"])
set𝝭!(elements["Γˡ"])
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
ops[5](elements["Γᵉ"],kᵇ,f)
ops[6](elements["Γˡ"],kᵇ,f)
ops[7](elements["Γᵇ"],kᵇ,f)

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
# e = abs(w[1]-𝑣)
# index = [8,16,32,64]
# XLSX.openxlsx("./xlsx/Circular_SSUL.xlsx", mode="rw") do xf
#     Sheet = xf[1]
#     ind = findfirst(n->n==ndivs,index)+1
#     Sheet["F"*string(ind)] = log10(5/ndiv)
#     Sheet["G"*string(ind)] = w
#     Sheet["H"*string(ind)] = log10(e)
# end

ind = 100

xs₁ = zeros(ind,ind)
ys₁ = zeros(ind,ind)
zs₁ = zeros(ind,ind)
xs₂ = zeros(ind,ind)
ys₂ = zeros(ind,ind)
zs₂ = zeros(ind,ind)
xs₃ = zeros(ind,ind)
ys₃ = zeros(ind,ind)
zs₃ = zeros(ind,ind)
xs₄ = zeros(ind,ind)
ys₄ = zeros(ind,ind)
zs₄ = zeros(ind,ind)
𝗠 = zeros(21)

for i in 1:ind
    for j in 1:ind
        Δy = (ind-i)*5/(ind-1)
        Δx = (ind-j)*(25-Δy^2)^0.5/(ind-1)
        x₁ = Δx
        y₁ = Δy
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

 for I in 1:ind
    for J in 1:ind
        xs₂[I,J] = -xs₁[ind-I+1,J]
        ys₂[I,J] = ys₁[ind-I+1,J]
        zs₂[I,J] = zs₁[ind-I+1,J]
        xs₃[I,J] = -xs₁[ind-I+1,ind-J+1]
        ys₃[I,J] = -ys₁[ind-I+1,ind-J+1]
        zs₃[I,J] = zs₁[ind-I+1,ind-J+1]
        xs₄[I,J] = xs₁[I,ind-J+1]
        ys₄[I,J] = -ys₁[I,ind-J+1]
        zs₄[I,J] =zs₁[I,ind-J+1]
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
surface!(xs₁,ys₁,zeros(ind,ind),color=zs₁,colorrange=(-2.5,0.13),shading = false,colormap=:rainbow)
surface!(xs₂,ys₂,zeros(ind,ind),color=zs₂,colorrange=(-2.5,0.13),shading = false,colormap=:rainbow)
surface!(xs₃,ys₃,zeros(ind,ind),color=zs₃,colorrange=(-2.5,0.13),shading = false,colormap=:rainbow)
surface!(xs₄,ys₄,zeros(ind,ind),color=zs₄,colorrange=(-2.5,0.13),shading = false,colormap=:rainbow)
# contour!(zs₁,levels=-2.5:0.1:0.13,color=:azure)
# Colorbar(fig[1,2], limits=(-2.5,0.13), colormap=:lightrainbow)


# save("./png/Circular_tri3_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/Circular_tri6_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
# save("./png/Circular_quad4_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)
save("./png/Circular_quad8_q1_"*string(ndiv)*"_"*string(ndivs)*".png",fig, px_per_unit = 3.0)

fig
