
import BenchmarkExample: BenchmarkExample
n = 8
# n₁ = 30
# n₂ = 29
# BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate/SquarePlate_"*string(n₁)*"_"*string(n₂)*".msh", transfinite = (n₁+1,n₂+1))
# BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate/SquarePlate_"*string(n)*".msh", transfinite = (n+1,n+1))
BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate/SquarePlate_tri6_"*string(n)*".msh", transfinite = n+1, order = 2)
# BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate/SquarePlate_quad_"*string(n)*".msh", transfinite = (n+1,n+1), quad = true)
# BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate/SquarePlate_quad_"*string(n₁)*"_"*string(n₂)*".msh", transfinite = (n₁+1,n₂+1), quad = true)
# BenchmarkExample.Circular.generateMsh("./msh/circular_"*string(n)*".msh", transfinite = n)
# BenchmarkExample.Circular.generateMsh("./msh/circular_quad_"*string(n)*".msh", transfinite = n, quad = true)
# BenchmarkExample.MorleysAcuteSkewPlate.generateMsh("./msh/MorleysAcuteSkewPlate_"*string(n)*".msh", transfinite = (n/2+1,n/2+1))1
# BenchmarkExample.MorleysAcuteSkewPlate.generateMsh("./msh/MorleysAcuteSkewPlate_"*string(n₁,n₂)*".msh", transfinite = (n₁/2+1,n₂/2+1))
