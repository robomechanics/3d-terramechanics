%Calculates normal stresses under wheel

function [sig, params] = sigma(params, theta)
%Unpack variables needed
theta_f = params.terr.theta_f;
theta_r = params.terr.theta_r;
theta_m = params.terr.theta_m;
k = params.soil.k;
n0 = params.soil.n0;
n1 = params.soil.n1;
n2 = params.soil.n2;

w = params.state.w;
v_x = params.state.v_x;
v_y = params.state.v_y;
V = sqrt(v_x^2 + v_y^2);
r_s = params.rover.r_s;

if w*r_s >= V %slip
    N = n0 + n1*(w*r_s - V)/(w*r_s);
else %skid
    N = n0 - n2*(w*r_s - V)/(V); %Note sign change on skid ratio to keep signs of all ns consistent
end

r = params.rover.r;
r_s = params.rover.r_s;
mu = params.rover.mu;

%Calculate normal stress at angle theta along the wheel
if (theta <= theta_f) && (theta >= theta_m) %Stresses in front of max stress
    sig_r = r^N*k*(cos(theta) - cos(theta_f))^N;
    sig_rs = r_s^N*k*(cos(theta) - cos(theta_f))^N;
elseif (theta < theta_m) && (theta >= theta_r) %Stresses behind max stress
    theta_e = theta_f - (theta - theta_r)/(theta_m - theta_r)*(theta_f - theta_m);
    sig_r = r^N*k*(cos(theta_e) - cos(theta_f))^N;
    sig_rs = r_s^N*k*(cos(theta_e) - cos(theta_f))^N;
else
    sig_r = 0;
    sig_rs = 0;
    fprintf("Why are you calculating stresses where the wheel isn't touching the soil dummy?\n")
end

sig = mu*sig_rs + (1-mu)*sig_r; %Average stress contributions from grousers and wheel surface
params.terr.sigma = sig;