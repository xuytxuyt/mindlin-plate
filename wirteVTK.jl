VTK_mix_pressure = quote
  #number of an elemnt nodes
  nₑₙ=3
  fo = open("./vtk/SquarePlate_mix_pressure_"*string(ndiv)*"_"*string(ndivs)*"Q2.vtk","w")
  @printf fo "# vtk DataFile Version 2.0\n"
  @printf fo "SquarePlate_mix\n"
  @printf fo "ASCII\n"
  @printf fo "DATASET POLYDATA\n"
  @printf fo "POINTS %i float\n" nˢ
  for p in nodes_s
    @printf fo "%f %f %f\n" p.x p.y p.z
  end
  @printf fo "POLYGONS %i %i\n" nₑₛ (nₑₙ+1)*nₑₛ
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
  @printf fo "POINT_DATA %i\n" nˢ
  @printf fo "SCALARS P float 1\n"
  @printf fo "LOOKUP_TABLE default\n"
  for p in nodes_s
    # @printf fo "%f\n" p.q₁ 
    @printf fo "%f\n" p.q₂
  end
  close(fo)
end