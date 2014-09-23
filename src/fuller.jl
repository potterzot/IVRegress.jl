
type FullerResp <: ModResp

  function FullerResp
  end
end

type FullerModel <: IVPredModel
  rr: FullerResp
  pp: IVPred
  fit::Bool
end



