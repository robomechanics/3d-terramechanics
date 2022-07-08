%Calculates tangential component of shear stresses under wheel
function [tau, params] = tau_t(params, theta)
c = params.soil.c;
sig = sigma(params, theta);
phi = params.soil.phi;
K = params.soil.K;
v_y = params.state.v_y;

vjl = -v_y;

vjt = v_jt(params,theta);
jt = j_t(params, theta);
jl = j_l(params, theta);

temp = (1 - exp(-sqrt(jt^2 + jl^2)/K));
if temp > 1
    temp = 1; %tau = tau_max
end

tau = sign(jt)*abs(vjt)/sqrt(vjt^2+vjl^2)*(c + sig*tan(phi))*temp;
if vjt == 0
    tau = 0; %case where theta=theta_0, otherwise we get NaN if vjy is also 0 (when beta = 0)
end
params.terr.tau_t = tau;