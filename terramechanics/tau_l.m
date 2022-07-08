%Calculates lateral component of shear stresses under wheel
function [tau, params] = tau_l(params, theta)
c = params.soil.c;
sig = sigma(params, theta);
phi = params.soil.phi;
K = params.soil.K;
v_y = params.state.v_y;

jt = j_t(params, theta);
jl = j_l(params, theta);

vjl = -v_y;

vjt = v_jt(params,theta);

temp = (1 - exp(-sqrt(jt^2 + jl^2)/K));
if temp > 1
    temp = 1; %tau = tau_max
end


tau = sign(jl)*abs(vjl)/sqrt(vjl^2+vjt^2)*(c + sig*tan(phi))*temp;
if vjl == 0
    tau = 0; %case where theta=theta_0 and beta=0 -> vjx = 0, which would give NaN
end

params.terr.tau_l = tau;