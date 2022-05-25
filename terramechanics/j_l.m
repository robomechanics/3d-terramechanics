%Compute lateral soil deformation module
function [j, params] = j_l(params, theta)

beta = params.state.beta;
b = params.rover.b;
theta_f = params.terr.theta_f;
theta_r = params.terr.theta_r;
w = params.state.w;
v_y = params.state.v_y;


v_jy = v_y;

j = (theta_f - theta)*v_jy/w; %classic defn
%to avoid 0*0/0 = Inf, manually assign those values:
if v_jy == 0
    j = 0;
elseif theta == theta_f
    j = 0;
end

params.terr.j_l = j;


