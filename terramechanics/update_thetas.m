%Gets entrance angle, exit angle, and maximum stress angle for the sinkage
%given by h and tuning parameters specified. 

function params = update_thetas(params)
r_s = params.rover.r_s;
h = params.terr.h;
v_x = params.state.v_x;
v_y = params.state.v_y;
V = sqrt(v_x^2 + v_y^2);
w = params.state.w;
a0 = params.tuned.a0;
a1 = params.tuned.a1;
b0 = params.tuned.b0;
b1 = params.tuned.b1;


theta_f = acos(1-h/r_s); 

if w*r_s >= V %slip
    theta_r = theta_f*(b0 + b1*(w*r_s - V)/(w*r_s));
    theta_m = theta_f*(a0 + a1*(w*r_s - V)/(w*r_s)); 
else %skid
    theta_r = params.tuned.theta_r0; 
    theta_m = params.tuned.theta_m0; 
end

if params.options.trench_on %Replace exit angle with the one found by trench model if option is on
    [~,~,~,~,~,~,~,~,ip] = trench_profile(params); % Use trench shape model to get exit angle, returns intersection point of soil profile and wheel centerline at the rear as 2d coordinates
    h_r = max(-ip(1,2), 0); %the sinkage at the exit angle of the wheel is defined as positive even though the coordinate is negative,
    theta_r = -acos(1+(h_r-h)/r_s);
end

c1 = theta_f/2; %TODO replace with a more physical approx
theta_0 = theta_f + (w*r_s - V)/(V)*(theta_f - c1); 

params.terr.theta_f = theta_f;
params.terr.theta_r = theta_r;
params.terr.theta_m = theta_m;
params.terr.theta_0 = theta_0;