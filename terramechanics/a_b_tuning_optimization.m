%Run a_b_tuning_optimization_wrapper and call this function to tune the
%values for a and b on the data. 
function [x, fval, exitflag] = a_b_tuning_optimization(all_results, h0, min_h_skid, x0, rovername)

% Optimization-based tuning of a0, a1, b0, and b1

%Set up parameters
params = initialize_params();
params = set_rover_params(params, rovername); %sets rover params 
params = set_soil_params(params, rovername); %sets soil and tuned params
params.options.trench_on = 0;
params.options.ishigami_sidewall = 0;

theta_f0 = acos(1 - h0/params.rover.r_s);
max_theta_m0 = acos(1 - min_h_skid/params.rover.r_s);
params.tuned.theta_f0 = theta_f0;


%Constraints and bounds
A = [1 1 0 0 0 0 0 0 0 0 0; 0 0 -1 -1 0 0 0 0 0 0 0];
b = [1 1]'; % Ax <= b -> a0 + a1 <=1, -b0 - b1 <= 1
lb = [0 0 -1 -1 0 -theta_f0 0 0 0 0 0]'; %Lower bounds
ub = [1 1 0 0 max_theta_m0 0 1.0 10^6 Inf Inf Inf]'; %Upper bounds
Aeq = [theta_f0 0 0 0 -1 0 0 0 0 0 0; 0 0 theta_f0 0 0 -1 0 0 0 0 0]; % 
beq = [0 0]'; % Ax = b -> a0*theta_f0 - theta_m = 0, b0*theta_f0 - theta_r = 0
nonlcon = []; %No nonlinear constraints

% if isempty(x0) %if we haven't passed in starting values, use these
%    x0 = [.2108 .6675 -.1207 -.2512 .1500 -.0859 .0135 2.3026 2.4215 .0924 2.2121]'; %Starting guesses, x = [a0 a1 b0 b1 theta_m0 theta_r0 K c n0 n1 n2]'
% end
%options = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'iter','OptimalityTolerance', 1e-16, 'FunctionTolerance', 1e-16, 'StepTolerance', 1e-16, 'DiffMinChange', 1e-16, 'FiniteDifferenceStepSize', 1e-6); %Saved as a comment as a record of the settings I have tried playing with


options = optimoptions(@fmincon, 'Algorithm', 'interior-point', 'Display', 'iter', 'ConstraintTolerance', 1e-6, 'OptimalityTolerance', 1e-3, 'FiniteDifferenceStepSize', 1e-9, 'UseParallel', true); %Sets the options for fmincon, algorithm choice shouldn't affect much
[x, fval, exitflag] = fmincon(@(x)a_b_tuning_objective_fn(x, params, all_results), x0, A, b, Aeq, beq, lb, ub, nonlcon, options);


