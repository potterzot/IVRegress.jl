#=
Implementation of LIML regressions with Bekker adjustments to standard errors

=#
module Kclass
using DataFrames

export 
  kclass #general k-class estimator interface






abstract IVModel <: RegressionModel



include("kclass.jl")





end # module
