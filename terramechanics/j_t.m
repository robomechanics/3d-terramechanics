%Compute tangential soil deformation module
function [j, params] = j_t(params, theta)
theta_f = params.terr.theta_f;
theta_r = params.terr.theta_r;
theta_0 = params.terr.theta_0;
v_x = params.state.v_x;
v_y = params.state.v_y;
V = sqrt(v_x^2 + v_y^2);
r_s = params.rover.r;
w = params.state.w;

if (r_s*w >= V) %slip condition
    j = r_s*(theta_f - theta) - (sin(theta_f) - sin(theta))*v_x/w;
else %skid condition
    if abs(theta - theta_0) < 1e-6
        j = 0; %This is always true in the above equations except for w=0, where we get NaN
    elseif (theta > theta_0) && (theta < theta_f) %front region
        j = v_x/w*((theta_f - theta)/(theta_f - theta_0)*(sin(theta_f) - sin(theta_0)) - (sin(theta_f) - sin(theta)));
    elseif (theta <= theta_0) && (theta > theta_r) %rear region
        j = r_s*(theta_0 - theta) - (sin(theta_0) - sin(theta))*v_x/w;
    else
        j = 0; %edge angles
    end
end
params.terr.j_t = j;

if isnan(j)
    fprintf('Check tolerances on theta = theta_0\n')
end
