
type TSLSResp <: ModResp

  function TSLSResp
  end
end

type TSLSModel <: IVPredModel
  rr: TSLSResp
  pp: IVPred
  fit::Bool
end



