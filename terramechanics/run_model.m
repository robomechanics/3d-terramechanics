%Takes in already averaged force and sinkage data, desired wheel and soil
%names, and runs the model with the corresponding values

function [all_params, Fx_approx, Fy_approx, h, sinkage_found] = run_model(wheel, betas, Vrys, trench_on, ishigami_sidewall, Fz_data, V_data, Vry_data, a0s, a1s, b0s, b1s)

%Initialize parameters with input options
params = initialize_params();
params = set_rover_params(params, wheel); %sets rover params
params = set_soil_params(params, wheel); %sets soil and tuned params
params.options.trench_on = trench_on;
params.options.ishigami_sidewall = ishigami_sidewall;
r = params.rover.r;

%Initialize output variables
Fx_approx = zeros(size(Fz_data));
Fy_approx = zeros(size(Fz_data));
Fz_approx = zeros(size(Fz_data));
sinkage_found = zeros(size(Fz_data));
h = zeros(size(Fz_data));
all_params = [];

%Calculate theta_m0, theta_r0 for each slip angle
theta_m0s = zeros(size(Fz_data));
theta_r0s = zeros(size(Fz_data));
for k=1:length(betas)
    V = V_data(1,k)/1000;
    params.state.beta = deg2rad(betas(k)); %We don't use the actual measured angle because they are very very close, but the case where we get beta sliiightly less than 0 breaks things
    beta = params.state.beta;
    params.state.w = V/params.rover.r_s;          %Set slip ratio to 0 [rad/s]
    params.state.v_x = V*cos(beta);               %Forward velocity of wheel [m/s]
    params.state.v_y = -V*sin(beta);               %Sideways velocity of wheel [m/s]

    %set tuned a and b values
    params.tuned.a0 = a0s(1,k);
    params.tuned.a1 = a1s(1,k);
    params.tuned.b0 = b0s(1,k);
    params.tuned.b1 = b1s(1,k);

    %set z force and use that to find sinkage and calculate other forces
    params.rover.W = Fz_data(4,k)/1000;             %Use the load applied in the test where s is approx 0, otherwise theta_m0 and theta_f0 will not be constant for a given slip ratio since the load varies slightly
    [~, params] = find_sinkage(params);
    params = update_thetas(params);

    %Set theta_m0 and theta_r0
    theta_m0s(:,k) = params.terr.theta_m;
    theta_r0s(:,k) = params.terr.theta_r;
end

%Compute model forces
for j=1:length(Vrys)
    for k=1:length(betas)
        %set state params
        V = V_data(j,k)/1000;
        params.state.beta = deg2rad(betas(k)); %We don't use the actual measured angle because they are very very close, but the case where we get beta sliiightly less than 0 breaks things
        beta = params.state.beta;
        params.state.w = Vry_data(j,k)/1000/r;       %Rotational (angular) velocity of wheel [rad/s]
        params.state.v_x = V*cos(beta);               %Forward velocity of wheel [m/s]
        params.state.v_y = -V*sin(beta);               %Sideways velocity of wheel [m/s]

        %set tuned a and b values
        params.tuned.a0 = a0s(j,k);
        params.tuned.a1 = a1s(j,k);
        params.tuned.b0 = b0s(j,k);
        params.tuned.b1 = b1s(j,k);
        params.tuned.theta_m0 = theta_m0s(j,k);
        params.tuned.theta_r0 = theta_r0s(j,k);

        %set z force and use that to find sinkage and calculate other forces
        params.rover.W = Fz_data(j,k)/1000;
        [h(j,k), params] = find_sinkage(params);
        sinkage_found(j,k) = params.terr.sinkage_found;
        params.terr.h = h(j,k);
        params = update_thetas(params);

        % compute forces on the wheel
        [Fx_approx(j,k), Fy_approx(j,k), Fz_approx(j,k), params] = forces(params);
        all_params = [all_params, params];
    end
end

end