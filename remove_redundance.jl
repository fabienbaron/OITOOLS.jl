using NearestNeighbors

function rm_redundance_slow(uv,uvtol)
  new_uv = uv[:,1];
  uv_length = size(uv,2);
  uv[:,find(uv[1,:].<0)] *= -1; # change any (-u,v) -> (u,-v)
  for i = 2:uv_length
    redundance = any(sum((broadcast(-,new_uv,uv[:,i])).^2,1).<uvtol^2);
    if (redundance == false)
      new_uv = hcat(new_uv,uv[:,i]);
    end
  end
  return new_uv
end

function rm_redundance_kdtree(uv,uvtol)
  @inbounds uv[:,find(uv[1,:].<0)] *= -1; # change any (-u,v) -> (u,-v)
  indx_redundance = collect(1:size(uv,2));
  kdtree = KDTree(uv);
  for (index,value) in enumerate(indx_redundance)
    redundance = inrange(kdtree,uv[:,value],uvtol);
    @inbounds indx_redundance[redundance] = minimum(redundance);
  end
  uv = uv[:,unique(indx_redundance)];
  return uv, indx_redundance;
end
