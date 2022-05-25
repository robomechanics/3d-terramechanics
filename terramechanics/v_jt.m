%Calculates tangential shear velocity at a point along the wheel surface
function v_jt = v_jt(params, theta)
v_x = params.state.v_x;
v_y = params.state.v_y;
V = sqrt(v_x^2 + v_y^2);
r_s = params.rover.r_s;
w = params.state.w;

theta_f = params.terr.theta_f;
theta_0 = params.terr.theta_0;
theta_r = params.terr.theta_r;

%determine shear velocity
if (r_s*w >= V) %slip condition
    v_jt = r_s*w - v_x*cos(theta);
else %skid condition
    if abs(theta - theta_0) < 1e-6
        v_jt = 0;
    elseif (theta > theta_0) && (theta <= theta_f) %front region
        v_jt = r_s*w - v_x*cos(theta) + ((sin(theta_f) - sin(theta_0))/(theta_f - theta_0)*v_x - r_s*w);
    elseif (theta <= theta_0) && (theta >= theta_r) %rear region
        v_jt = r_s*w - v_x*cos(theta);
    else
        v_jt = 0; %edge angles, not in contact
    end
end
