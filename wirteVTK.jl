VTK_mix_pressure = quote
    #number of an elemnt nodes
    nₑₙ=3
    fo = open("./vtk/SquarePlate_mix_pressure_"*string(ndiv)*"_"*string(ndivs)*"Q1.vtk","w")
    @printf fo "# vtk DataFile Version 2.0\n"
    @printf fo "SquarePlate_mix\n"
    @printf fo "ASCII\n"
    @printf fo "DATASET POLYDATA\n"
    @printf fo "POINTS %i float\n" nᵇ
    
  for p in nodes
    @printf fo "%f %f %f\n" p.x p.y p.z
  end
    @printf fo "POLYGONS %i %i\n" nₑ (nₑₙ+1)*nₑ
    for ap in elements["Ω"]
       𝓒 = ap.𝓒
       if nₑₙ==3
        @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==4
        @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==6
        @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       elseif nₑₙ==8
        @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
       end
    end
    # @printf fo "POINT_DATA %i\n" nᵤ
    @printf fo "CELL_DATA %i\n" nₑ
    @printf fo "SCALARS PRESSURE float 1\n"
    @printf fo "LOOKUP_TABLE default\n"
      for ap in elements["Ω"]
        𝓒 = ap.𝓒
        𝓖 = ap.𝓖
        γ₁ = 0.0 
        γ₂ = 0.0
        for (i,ξ) in enumerate(𝓖)
                N = ξ[:𝝭]
                B₁ = ξ[:∂𝝭∂x]
                B₂ = ξ[:∂𝝭∂y]
                for (j,xⱼ) in enumerate(𝓒)
                    γ₁ += B₁[j]*xⱼ.d₁ - N[j]*xⱼ.d₂
                    γ₂ += B₂[j]*xⱼ.d₁ - N[j]*xⱼ.d₃
                end
               
        end
        p = γ₁
        @printf fo "%f\n" p   
    end
    close(fo)
    end

    VTK_mix_pressure_E = quote
      #number of an elemnt nodes
      nₑₙ=3
      fo = open("./vtk/SquarePlate_mix_pressure_"*string(ndiv)*"_"*string(ndivs)*"EQ1.vtk","w")
      @printf fo "# vtk DataFile Version 2.0\n"
     @printf fo "SquarePlate_mix\n"
      @printf fo "ASCII\n"
      @printf fo "DATASET POLYDATA\n"
      @printf fo "POINTS %i float\n" nₚ
     
      for p in nodes_p
         @printf fo "%f %f %f\n" p.x p.y p.z
      end
      @printf fo "POLYGONS %i %i\n" nₑₚ (nₑₙ+1)*nₑₚ
      for ap in Ω
         𝓒 = ap.𝓒
         if nₑₙ==3
          @printf fo "%i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==4
          @printf fo "%i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==6
          @printf fo "%i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         elseif nₑₙ==8
          @printf fo "%i %i %i %i %i %i %i %i %i\n" nₑₙ (x.𝐼-1 for x in 𝓒)...
         end
      end
      @printf fo "POINT_DATA %i\n" nₚ
      @printf fo "SCALARS P float 1\n"
      @printf fo "LOOKUP_TABLE default\n"
      for p in nodes_p
         @printf fo "%f\n" p.Q₁ 
      end