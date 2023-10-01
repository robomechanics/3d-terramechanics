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
    % Hand tuned for grousered wheel, goes with trench_off_grouser.mat
    params.soil.K = 0.021;                  %Shear modulus [m]
    params.soil.c = 1.00;                   %Cohesion stress [kPa]
    params.soil.n0 = 1.46;                  %Sinkage exponent constant []
    params.soil.n1 = 0.01;                 %Sinkage exponent constant []
    params.soil.n2 = .74;                   %Sinkage exponent constant []
    %Linear tuning variables for the grousered wheel, fitted from tuned
    %results
    params.tuned.a_c = 0.2365;
    params.tuned.a_b = 0.007744;
    params.tuned.a_s = 0.6892;
    params.tuned.a_sb = -0.007491;
    params.tuned.b_c = [];
    params.tuned.b_b = [];
    params.tuned.b_s = [];
    params.tuned.b_sb = [];
    fprintf('Using hand tuned grousered wheel soil values\n')
elseif strcmp(soil, 'K10_mini_smooth')
    % Hand tuned for smooth wheel, goes with trench_off_smooth.mat
    params.soil.K = 0.021;                   %Shear modulus [m]
    params.soil.c = 1.0;                   %Cohesion stress [kPa]
    params.soil.n0 = 1.46;                   %Sinkage exponent constant []
    params.soil.n1 = .01;                   %Sinkage exponent constant []
    params.soil.n2 = .55;                   %Sinkage exponent constant []
    %Linear tuning variables for the smooth wheel, fitted from tuned
    %results
    params.tuned.a_c = 0.4045;
    params.tuned.a_b = 0.005841;
    params.tuned.a_s = 0.4198;
    params.tuned.a_sb = -0.004566;
    params.tuned.b_c = [];
    params.tuned.b_b = [];
    params.tuned.b_s = [];
    params.tuned.b_sb = [];
    fprintf('Using hand tuned smooth wheel soil values with r_s=1.02r\n');
elseif strcmp(soil, 'K10_mini_legacy_ICRA')
    % Old values used in ICRA 2019 paper
    params.soil.K = .025;                   %Shear modulus [m]
    params.soil.c = 0;                   %Cohesion stress [kPa]
    params.soil.n0 = 1;                   %Sinkage exponent constant []
    params.soil.n1 = 0;                  %Sinkage exponent constant []
    params.soil.n2 = 0;                   %Sinkage exponent constant []
    params.soil.phi = deg2rad(29);        %Angle of repose [rad]
    params.tuned.a0 = .4;
    params.tuned.a1 = .15;
    params.tuned.b0 = 0;
    params.tuned.b1 = 0;
    params.soil.k = 1523.4 + 0.9/.048; % Sinkage modulus, specific to wheel width of .048m
    fprintf('Using legacy soil values from ICRA 2019 paper');
elseif strcmp(soil, 'KRex_roverscape')
    % Soil in NASA Ames Roverscape fall 2022: crushed granite, poorly
    % graded (ranges from dust to >1cm pebbles, majority of grains are >2mm
    params.soil.k = 1203.54;                        %Sinkage modulus [kN/m^(n+m+2)] value from Tsubaki & Ishigami 2021
    params.soil.k = 2500;
    params.soil.rho = 1705*9.8/10^3;       %Soil specific weight [kN/m^3] from density, measured from loose soil
    %params.soil.rho = 1963*9.8/10^3;       %Soil specific weight [kN/m^3] from density, measured from manually tamped soil
    params.soil.phi = deg2rad(38);              %Angle of repose [rad], measured from loose piles at several points in Roverscape
    params.soil.K = .04;                    %Shear modulus [m] guesstimated from typical range of 10-25 cm per Wong TGV 2022
    params.soil.c = 5.0;                  %Cohesion stress [kPa] value from Tsubaki & Ishigami 2021
    params.soil.n0 = 1.7;                   %Sinkage exponent constant [] value from Tsubaki & Ishigami 2021
    params.soil.n0 = .6;
    params.soil.n1 = .5;    %guess              %Sinkage exponent constant []
    params.soil.n2 = 0.7;    %guess               %Sinkage exponent constant []
    params.tuned.a_c = 0.40;                  %copied from smooth wheel above
    params.tuned.a_b = 0.006;                  %copied from smooth wheel above
    params.tuned.a_s = 0.42;                  %copied from smooth wheel above
    params.tuned.a_sb = -0.005;                  %copied from smooth wheel above
    params.tuned.b_c = [];                  %copied from smooth wheel above
    params.tuned.b_b = [];                  %copied from smooth wheel above
    params.tuned.b_s = [];                  %copied from smooth wheel above
    params.tuned.b_sb = [];                  %copied from smooth wheel above
    fprintf('Using hand tuned soil values for KREX2 in the roverscape with crushed granite substrate\n');
else
    error('Error - valid soil name not specified. Use either "K10_mini_grouser" or "K10_mini_smooth". \n');
end

%Destructive angle has to be defined after phi!
params.soil.X_c = 45*pi/180 - params.soil.phi/2; %destructive angle formulation for when surcharge is fully formed, Wong 2008 [rad]

vars = params;