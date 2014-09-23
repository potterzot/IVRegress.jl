using StatsBase
using Base.Test

# Run the tests
function test()
  #=
  Test the LIML implementation using mroz dataset. Results should be:

  lwage is dependent
  hours is instrumented
  exper is non-instrumented
  educ, kidslt6, kidsge6, age, nwifeinc are instruments (include also any non-instrumented)
 
  1st-Stage


 
  2nd-Stage
          Coef      Std Err   Bekker  'Small'
  hours    .001789  .001115   .00203  .001148
  exper   -.092979  .141521   .25486  .145624
  =#
    
  mroz = readtable("/home/potterzot/data/stata/mroz.csv") 
 
  # Drop incomplete cases
  data = complete_cases!(mroz)

  # Small sample for testing
  data = data[data[:age].<=31, :]

  
  # Dependent variable 
  Y = data[:lwage]

  # Instrumented Variable
  X = data[:hours]

  # Instruments
  Z = hcat(data[:educ], data[:kidslt6], data[:kidsge6], data[:age], data[:nwifeinc])

  result = iv(Y, X, Z)
end

