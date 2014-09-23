#=
Two Stage Least Squares (2SLS)
Limited Information Maximum Likelihood (LIML)
Generalized Method of Moments (GMM)

Required Variables:
y is dependent variable
X is instrumented variables (endogenous variables)
Z is instruments (excluded exogenous)
W is non-instrumented exogenous variables

NOTE: All "W" non-instrumented variables should, by definition, be included in the first-stage regression. This is the behavior by default.

Options:
wt is an NxN matrix of weights
bekker is a bekker adjustment to standard errors

In stata form:

ivreg2 liml Y W (X=Z)

In equation form:
Y = X*b + W*d + u
X = Z*a + v

Formula is:
y ~ W + (X = 

=#
type LIMLModel <: IVModel
  rr: LIMLResp
  pp: IVPred
  fit::Bool
end
type FullerModel <: IVModel
  rr: FullerResp
  pp: IVPred
  fit::Bool
end



type KclassResp{V<:FPVector} <: ModResp
  y::V
  var::V
  wts::V
  alpha::float
  fuller::float
  function KclassResp() #constructor
    res = new() #something goes here...
    updateemu!(res, eta)
  end
end

function kclass(f::Formula, df::AbstractDataFrame, model::String="2SLS", fuller::Int=0)
  using StatsBase: StatsBase, CoefTable, StatisticalModel, RegressionModel
  import StatsBase: coef, coeftable, confint, deviance, loglikelihood, nobs, stderr, vcov, residuals, predict, fit
  
  mf = ModelFrame(f, df) 
  mm = ModelMatrix(mf)
  rr = 
  res = KclassMod(mf, rr, 


  vcov(x::KclassModel) = scale(x, true) * inv(cholfact(x.pp))
  stderr(x::KclassModel) = sqrt(diag(vcov(x))) 

  # pinv is equivalent but slower to inv(cholfact(Z'*Z))


  # Number of each type of variable
  N = size(y, 1) #number of obs
  K = size(X, 2) # 
  L = size(W, 2)
  M = size(Z, 2)

  df_ess = K + L
  df_rss = N - df_ess

  #Calculate the matrixMatrices
  A = hcat(y,X,W,Z)
  Y = hcat(y,X)
  AA = A' * A
  
  #Sub matrices
  yy = AA[1,1] # 1x1
  XX = AA[2:1+K, 2:1+K] # KxK
  Xy = AA[2:1+K, 1] # Kx1
  WW = AA[2+K:2+K+L, 2+K:2+K+L] # LxL
  ZZ = AA[2+K+L:end, 2+K+L:end] # MxM
  XZ = AA[2:1+K, 2+K:end] # KxM
  PZ = Z * pinv((Z' * Z) * Z' # NxN
  MZ = eye(N) - PZ # NxN
  XPZX = X' * PZ * X # KxK
  XPZy = X' * PZ * y # Kx1 
  YMY = Y' * MZ * Y # 1+K x 1+K
  YPZY = YY - YMY

  #Get lambda/alpha
  if model=="2sls"
    alpha = 0
  elseif(model=="liml" || model=="fuller")
    HH = YY * pinv(YPZY)
    lambda = eigmin(HH) # lambda is 1.84247. alpha = 1 - 1/lambda = .4572
    if model=="liml"
      alpha = 1 - 1/lambda
    elseif model=="fuller"
      alpha = ((1 - 2/lambda) * fuller / N) / ((1 + 1/lambda) * fuller / N) #C is fuller param 
  end
  ###BETAS
  # first regression (Z)
  betaZ = pinv(ZZ) * XZ'
 
  # second regression (X)
  H = XPZX - alpha * XX
  betaX = pinv(H) * (XPZy - alpha * Xy)
  
  # Variances
  err = Y - X * betaX #Residuals
  rss = err' * err #Sum of Squared Residuals / Residual SS
  tss = (y - mean(y))' * (y - mean(y)) #Total sum of squares
  mss = rss - tss # model sum of squares
  mse = rss / df_rss # Variance of Residuals, also sigma squared of the regression
  rmse = sqrt(mse) # standard error of the regression
  tss_uncentered = yy[1]
  r_squared_centered = 1-rss/tss # centered R-Squared
  r_squared_uncentered = 1-rss/tss_uncentered # uncentered R-squared
 
  #Z beta Variance
  err1 = X .- Z * betaZ
  mse1 = (err1' * err1) / (N-M)
  rmse1 = sqrt(mse1)
  seZ = rmse1 .* sqrt(diag(pinv(zz)))
 
  #X beta Variance
  Xbar = X - err .* (err' * X)/rss
  Sigma = mse * ( (1-alpha)^2 * (Xbar' * PZ * Xbar) + alpha^2 * (Xbar' * MZ * Xbar) )
  seX = rmse .* sqrt(diag(pinv(xx)))
  seX_bekker = sqrt(pinv(H) * Sigma * pinv(H))
  #T-statistic
  betaX / seX

  
  
  return betaX
end

#Variants on inputs
kclass(e::Expr) = kclass(Formula(e))
kclass(f::Formula) = kclass(f)
kclass(s::String) = kclass(Formula(parse(s)[1]))

function predict(mm::KclassMod, newx::Array)
  newx * coef(mm)
end


