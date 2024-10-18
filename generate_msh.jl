
import BenchmarkExample: BenchmarkExample
# n = 
n₁ = 9
n₂ = 7
BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate_"*string(n₁)*"_"*string(n₂)*".msh", transfinite = (n₁+1,n₂+1))
# BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate_"*string(n)*".msh", transfinite = (n+1,n+1))
# BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate_quad_"*string(n)*".msh", transfinite = n, quad = true)
# BenchmarkExample.Circular.generateMsh("./msh/circular_"*string(n)*".msh", transfinite = n)
# BenchmarkExample.Circular.generateMsh("./msh/circular_quad_"*string(n)*".msh", transfinite = n, quad = true)
# BenchmarkExample.MorleysAcuteSkewPlate.generateMsh("./msh/MorleysAcuteSkewPlate_"*string(n)*".msh", transfinite = (n/2+1,n/2+1))1
# BenchmarkExample.MorleysAcuteSkewPlate.generateMsh("./msh/MorleysAcuteSkewPlate_"*string(n₁,n₂)*".msh", transfinite = (n₁/2+1,n₂/2+1))