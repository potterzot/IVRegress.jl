type KclassResp <: ModResp

  function KclassResp
  end
end

type KclassModel <: IVPredModel
  rr: KclassResp
  pp: IVPred
  fit::Bool
end



