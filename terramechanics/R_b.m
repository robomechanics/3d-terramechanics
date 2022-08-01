%Forces on side face of the wheel 
%calculates the incremental bulldozing resistance for a portion of the
%side face of the wheel as a function of theta (goes into integral for 
%calculating side forces)

function [Rb, params] = R_b(params, theta)

%Get needed parameters
c = params.soil.c;
rho = params.soil.rho; %density of soil
r_s = params.rover.r_s;
phi = params.soil.phi; %friction angle, determines the height of soil piled in fornt of wheel
X_c = params.soil.X_c; %destructive angle, determines angle of sheared soil below surface
beta = params.state.beta;
theta_f = params.terr.theta_f;
theta_r = params.terr.theta_r;

hb = min(2*r_s*cos(theta), (sin(theta) - sin(theta_r))/(sin(theta_f) - sin(theta_r))*r_s*(cos(theta_r) - cos(theta_f)) - r_s*(cos(theta_r) - cos(theta))); 

if params.options.ishigami_sidewall
    hh = r_s*(cos(theta) - cos(theta_f));
    D1 = cot(X_c) + tan(X_c + phi);
    D2 = cot(X_c) + (cot(X_c)^2)/cot(phi);
    Rb = D1*(c*hh + D2*rho*hh^2/2);
else
    Rb = 1/2*rho*hb.^2*(cot(X_c)^2)*(1 + cot(X_c)*tan(phi)/2) + 2*hb.*c*cot(X_c); %Wong 2008
    Rb = Rb*abs(sin(beta)); 
end

%Don't want negative forces at high sinkage
if (theta>pi/2) || (theta < -pi/2)
    Rb = 0;
end
%Catch edge cases, since numerical error will result in hh !=0 at the edge
%of the contact angles when it should be 0
if (theta == params.terr.theta_f) || (theta == params.terr.theta_r) 
    Rb = 0;
end
Rb = Rb;
