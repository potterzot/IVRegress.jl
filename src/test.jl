function ivreg(form::Formula, data)#form::Formula,d) #df::Array)
#  using StatsBase: StatsBase, CoefTable, StatisticalModel, RegressionModel
#  import StatsBase: coef, coeftable, confint, deviance, loglikelihood, nobs, stderr, vcov, residuals, predict, fit

  #mf = ModelFrame(form, df) 
  #mm = ModelMatrix(mf)
  #return mm
#  return true
end  
#ivreg(form::Formula, df::Array) = ivreg(form, convert(DataFrame, df))

function test(a::Int, b::Int)
  return a*b
end

ivreg(form, df) = 1





