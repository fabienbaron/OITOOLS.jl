using SpecialFunctions
#uniform disk
#param[1] = diameter in mas
#v_r = radius
function visibility_ud(param, v_r;tol=1e-14)
theta = param[1]/2.0626480624709636e8;
V = 2.0*(besselj1.(pi*theta*v_r))./(pi*theta*v_r)
indx= findall(abs.(theta*v_r).<tol)
if indx !=[]
    V[indx].=1.0;
end
return V
end

#power law
function visibility_ldpow(param, v_r)
maxk=200;
theta = param[1]/2.0626480624709636e8;
#f=k->gamma(param[2]/2+2)/(gamma(param[2]/2+k+2)*gamma(k+1))*(-0.25*(pi*theta*v_r).^2).^k; # note: can lead to NaN due to typical gamma blow up
f=k->(-1)^k*exp.(lgamma(param[2]/2+2)-lgamma(param[2]/2+k+2)-lgamma(k+1).+2k*log.(0.5*pi*theta*v_r));
V = sum(map(f,collect(0:maxk)));
return V
end

#quadratic law
function visibility_ldquad(param,v_r)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*v_r*theta);
V = ((1.0-param[2]-param[3])*(besselj1.(zeta)./zeta)+((theta+2*param[3])/sqrt(2/pi))*(
(sqrt.(2 ./(pi*zeta)).*((sin.(zeta)./zeta)-cos.(zeta)))./zeta.^(3/2))-2*param[3]*
(besselj.(2, zeta)./zeta.^2))./(0.5-param[2]/6-param[3]/12)
return V
end

#linear law
function visibility_ldlin(param,v_r)
theta = param[1]/2.0626480624709636e8;
zeta = (pi*v_r*theta);
V = ((1.0-param[2])*(besselj1.(zeta)./zeta)+theta/sqrt(2/pi)*(
(sqrt.(2 ./(pi*zeta)).*((sin.(zeta)./zeta)-cos.(zeta)))./zeta.^(3/2)))./(0.5-param[2]/6)
return V
end

# visibility of an annulus of unit flux
function visibility_annulus(r_in, r_out, v_r;tol=1e-10)
if abs(r_out - r_in)<tol
  return zeros(Complex{Float64},length(v_r));
else
  return (1.0+0*im).*(visibility_ud([2*r_out],v_r)*r_out^2-visibility_ud([2*r_in],v_r)*r_in^2)/(r_out^2-r_in^2);
end
end
