
import BenchmarkExample: BenchmarkExample
n = 64
BenchmarkExample.SquarePlate.generateMsh("./msh/SquarePlate_"*string(n)*".msh", transfinite = n)
# BenchmarkExample.Circular.generateMsh("./msh/circular_"*string(n)*".msh", transfinite = n)
# BenchmarkExample.Circular.generateMsh("./msh/circular_quad_"*string(n)*".msh", transfinite = n, quad = true)
# BenchmarkExample.MorleysAcuteSkewPlate.generateMsh("./msh/MorleysAcuteSkewPlate_"*string(n)*".msh", transfinite = n)