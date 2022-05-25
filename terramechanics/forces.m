%Computes the forces on a wheel by numerically integrating stresses over
%the surface.
function [Fx, Fy, Fz, params] = forces(params)

%Soil and platform constants, these never change 
b = params.rover.b; %width of wheel in m
r_s = params.rover.r_s; %shearing radius of wheel (usually wheel radius + grouser height)
beta = params.state.beta; %angle of wheel relative to direction of travel

%Find sweep angles, used in integration limits. These are highly state dependent
params = update_thetas(params);
theta_range = linspace(params.terr.theta_r, params.terr.theta_f, 101);

%compute stresses and forces over wheel surface
sigmas = zeros(length(theta_range),1);
tau_ls = sigmas;
tau_ts = sigmas;
R_bs = sigmas;
hh = sigmas;
Fx_approx = 0;
Fy_approx = 0;
Fz_approx = 0;

for i=1:length(theta_range)
    theta = theta_range(i);
    [sigmas(i), params] = sigma(params, theta); %compute normal stress
    [tau_ls(i), params] = tau_l(params, theta); %compute lateral shear stress
    [tau_ts(i), params] = tau_t(params, theta); %compute tangential shear stress
    [R_bs(i), params] = R_b(params, theta); %compute bulldozing force on side face
    hh(i) = r_s*(cos(theta) - cos(params.terr.theta_f));
    if i>1     
        if params.options.ishigami_sidewall == 1
            Px = r_s*b*(-sigmas(i)*sin(theta) + tau_ts(i)*cos(theta));
            Fx_approx = Fx_approx + (theta-theta_range(i-1))*Px;
            Py = r_s*b*(tau_ls(i)) + (r_s - hh(i)*cos(theta))*R_bs(i)*sin(beta);
            Fy_approx = Fy_approx + (theta-theta_range(i-1))*Py;
        else
            Px = r_s*b*(-sigmas(i)*sin(theta) + tau_ts(i)*cos(theta));
            Fx_approx = Fx_approx + (theta-theta_range(i-1))*Px;
            Py = r_s*b*(tau_ls(i)) + r_s*cos(theta)*R_bs(i)*sin(beta);
            Fy_approx = Fy_approx + (theta-theta_range(i-1))*Py;
        end
        Pz = r_s*b*(sigmas(i)*cos(theta) + tau_ts(i)*sin(theta));
        Fz_approx = Fz_approx + (theta-theta_range(i-1))*Pz;
    end
end

    if any(isnan([Fx_approx, Fy_approx, Fz_approx]))||Fy_approx>0 ||any([imag(Fx_approx) ~= 0, imag(Fy_approx) ~= 0, imag(Fz_approx) ~= 0])
        Fx_approx;
    end

Fx = Fx_approx;
Fy = Fy_approx;
Fz = Fz_approx;