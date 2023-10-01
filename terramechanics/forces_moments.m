%Computes the forces on a wheel by numerically integrating stresses over
%the surface. Currently, the moment calculations neglect the contribution by
% cohesion in Rb.
function [Fx, Fy, Fz, Mx, My, Mz, params] = forces_moments(params)

%Soil and platform constants, these never change 
b = params.rover.b; %width of wheel in m
r_s = params.rover.r_s; %shearing radius of wheel (usually wheel radius + grouser height)
beta = params.state.beta; %angle of wheel relative to direction of travel
rho = params.soil.rho;
c = params.soil.c;
X_c = params.soil.X_c;
phi = params.soil.phi;


%Find sweep angles, used in integration limits. These are highly state dependent
params = update_thetas(params);
theta_range = linspace(params.terr.theta_r, params.terr.theta_f, 101);

%compute stresses and forces over wheel surface
sigmas = zeros(length(theta_range),1);
tau_ls = sigmas;
tau_ts = sigmas;
R_bs = sigmas;
hh = sigmas;
hb = sigmas;
h_COP = sigmas;
Fx_approx = 0;
Fy_approx = 0;
Fz_approx = 0;
Mx_approx = 0;
My_approx = 0;
Mz_approx = 0;

for i=1:length(theta_range)
    theta = theta_range(i);
    [sigmas(i), params] = sigma(params, theta); %compute normal stress
    [tau_ls(i), params] = tau_l(params, theta); %compute lateral shear stress
    [tau_ts(i), params] = tau_t(params, theta); %compute tangential shear stress
    [R_bs(i), params] = R_b(params, theta); %compute bulldozing force on side face
    hh(i) = r_s*(cos(theta) - cos(params.terr.theta_f)); %depth used in Ishigami calculation
    hb(i) = min(2*r_s*cos(theta), (sin(theta) - sin(params.terr.theta_r))/(sin(params.terr.theta_f) - sin(params.terr.theta_r))*r_s*(cos(params.terr.theta_r) - cos(params.terr.theta_f)) - r_s*(cos(params.terr.theta_r) - cos(theta))); %Bulldozing depth calc 
    h_COP(i) = (1/3*hb(i)^2*rho*cot(X_c)^2*(1+.5*cot(X_c)*tan(phi)) + 1/2*hb(i)*c*cot(X_c))/(1/2*hb(i)*rho*cot(X_c)^2*(1+.5*cot(X_c)*tan(phi)) + c*cot(X_c)); %Center of pressure, calculated using my formulation for bulldozing stress only
    if hb(i) == 0 && c == 0 %Catch 0/0 error at end points from 0 cohesion
        h_COP(i) = 0;
    end
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

        Mx_approx = Mx_approx - (r_s^2*b*tau_ls(i)*cos(theta) + R_bs(i)*sin(beta)*(r_s*cos(theta))*(h_COP(i)))*(theta-theta_range(i-1)); %TODO check this
        My_approx = My_approx - r_s^2*b*tau_ts(i)*(theta-theta_range(i-1)); %TODO check this
        Mz_approx = Mz_approx - (r_s^2*b*tau_ls(i)*sin(theta) + R_bs(i)*sin(beta)*r_s^2*sin(theta)*cos(theta))*(theta-theta_range(i-1)); %TODO check this
    end
end

if any(isnan([Fx_approx, Fy_approx, Fz_approx, Mx_approx, My_approx, Mz_approx])) ||any([imag(Fx_approx) ~= 0, imag(Fy_approx) ~= 0, imag(Fz_approx) ~= 0])
    Fx_approx;
end

Fx = Fx_approx;
Fy = Fy_approx;
Fz = Fz_approx;
Mx = Mx_approx;
My = My_approx;
Mz = Mz_approx;

params.terr.Fx = Fx;
params.terr.Fy = Fy;
params.terr.Fz = Fz;
params.terr.Mx = Mx;
params.terr.My = My;
params.terr.Mz = Mz;