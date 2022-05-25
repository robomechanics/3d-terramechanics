% Set the soil values needed to run terramechanics simulations. 
% User must specify soil name from the following list: 'K10_mini_grouser',
%'K10_mini_smooth'

function vars = set_soil_params(params, soil)
if ~exist('soil')
    error('Error - valid soil name not specified. Use either "K10_mini_grouser" or "K10_mini_smooth".');
end

params.soil.k = 8e3;                        %Sinkage modulus [kN/m^(n+m+2)] FROM TESTED VALUES (rounded from median grousered wheel value, rounded from mean and median plate values)
params.soil.rho = 1.33*10^3*9.8/10^3;       %Soil specific weight [kN/m^3] from density
params.soil.phi = deg2rad(29);              %Angle of repose [rad]
% params.soil.c = .0;                       %Cohesion stress [kPa]
% params.soil.K = .025;                     %Shear modulus [m]
% params.soil.n0 = params.soil.n;           %Sinkage exponent constant []
% params.soil.n1 = 0;                       %Sinkage exponent constant []
% params.soil.n2 = params.soil.n;           %Sinkage exponent constant []

if strcmp(soil, 'K10_mini_grouser')
    % Hand tuned for grousered wheel, goes with trench_off_grouser.mat, which has
    % exit flag 2 for beta = 15, 75 exit flag 1 for beta = 15, 30, 45, 60, 90
    % and trench_on_grouser.mat, which has
    % exit flag 2 for beta = 30, 60, 75, and exit flag 1 for beta = 0, 15, 45, 90
    params.soil.K = 0.021;                  %Shear modulus [m]
    params.soil.c = 1.14;                   %Cohesion stress [kPa]
    params.soil.n0 = 1.45;                  %Sinkage exponent constant []
    params.soil.n1 = 0.005;                 %Sinkage exponent constant []
    params.soil.n2 = .75;                   %Sinkage exponent constant []
    fprintf('Using hand tuned grousered wheel soil values\n')
elseif strcmp(soil, 'K10_mini_smooth')
    % Hand tuned for smooth wheel, goes with trench_off_smooth.mat, which has
    % exit flag 2 for beta = 30, 45, 60,  exit flag 1 for beta = 0, 15, 75, 90
    % and trench_on_smooth.mat, which has
    % exit flag 2 for beta = ~, exit flag 1 for beta = ~
    params.soil.K = 0.01;                   %Shear modulus [m]
    params.soil.c = 1.14;                   %Cohesion stress [kPa]
    params.soil.n0 = 1.4;                   %Sinkage exponent constant []
    params. soil.n1 = .01;                  %Sinkage exponent constant []
    params.soil.n2 = .55;                   %Sinkage exponent constant []
    fprintf('Using hand tuned smooth wheel soil values with r_s=1.02r\n');
else
    error('Error - valid soil name not specified. Use either "K10_mini_grouser" or "K10_mini_smooth". \n');
end

%Destructive angle has to be defined after phi!
params.soil.X_c = 45*pi/180 - params.soil.phi/2; %destructive angle formulation for when surcharge is fully formed, Wong 2008 [rad]

vars = params;