function error = a_b_tuning_objective_fn(x, params, all_results)
%all_results is the data to be tuned over, if you want to tune a and b for
%a subset of data you must crop it down before passing it to this function

%perform the params setup before passing it to this function - make sure to
%set rover and soil properties, but then nullify the values for a0, a1, b0
%and b1

%initialize the error to 0
error = 0;

%Set a and b values
params.tuned.a0 = x(1);
params.tuned.a1 = x(2);
params.tuned.b0 = x(3);
params.tuned.b1 = x(4);
params.tuned.theta_m0 = x(5);
params.tuned.theta_r0 = x(6);
% params.soil.K = x(7);
% params.soil.c = x(8);
% params.soil.n0 = x(9);
% params.soil.n1 = x(10);
% params.soil.n2 = x(11);

for i=1:length(all_results)
    %Extract measured forces to compare model prediction to
    datapt = all_results(i);
    avg_Fx = -datapt.avg_Fy/1000;
    avg_Fy = datapt.avg_Fx/1000;
    avg_Fz = -datapt.avg_Fz/1000;

    %Extract measured speeds, angles and sinkages to feed into the model
    %prediction
    V = datapt.avg_Vx/1000;
    w = datapt.avg_Vry/1000/params.rover.r;
    h = -double(datapt.avg_Z/1000);
    beta = deg2rad(datapt.beta); %Note that this uses the commanded angle rather than the measured angle, which is avg_angle


    %Set values for model to run on
    params.state.v_x = V*cos(beta);
    params.state.v_y = -V*sin(beta);
    params.state.w = w;
    params.state.beta = beta;
    params.terr.h = h;

    %Calculate angles
    params = update_thetas(params);

    %Compute forces
    [Fx, Fy, Fz, params] = forces(params);

    %Calculate value of objective function
    err = 1*(Fx - avg_Fx)^2/(avg_Fx^2) + 1*(Fy - avg_Fy)^2/(avg_Fy^2) + 1*(Fz - avg_Fz)^2/(avg_Fz^2);
    error = error + err;

    %         if imag(err) ~= 0
    %             %err = (imag(err)*10^10)
    %             err;
    %         end
end
error = double(error);



