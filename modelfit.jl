
function fit_model_v2(data, visfunc, init_param)
  nparams=length(init_param)
  indx= data.indx_v2
  nv2 = length(data.v2[indx]);
  r=data.v2_baseline
  chisq=(param,g)->sum(((abs2(visfunc(param,r)-data.v2[indx])./data.v2_err[indx]).^2))/nv2;
  opt = Opt(:LN_NELDERMEAD, nparams);
  min_objective!(opt, chisq)
  (minf,minx,ret) = optimize(opt, init_param);
  println("got $minf at $minx (returned $ret)")
  cvis_model = visfunc(minx,data.v2_baseline)
  return (minf,minx,ret,cvis_model)
end
