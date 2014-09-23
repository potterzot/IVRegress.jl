
type LIMLResp <: ModResp

  function LIMLResp
  end
end

type LIMLModel <: IVPredModel
  rr: LIMLResp
  pp: IVPred
  fit::Bool
end



