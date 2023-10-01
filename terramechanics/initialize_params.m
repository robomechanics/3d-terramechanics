%Initializes struct of structs with all values uswed in running the
%terramechanics model.

function params = initialize_params()

% set rover geometry (these values do not change throughout terramechanics
% evaluations)
r = [];                 %Wheel radius, measured to base of grousers [m]
b = [];                 %Wheel width [m]
h_g = [];               %Grouser height [m]
r_s = [];               %Shearing radius of soil [m]
W = [];                 %Load on one wheel [kN]
mu = [];                %Area ratio of grousers []
zeta = [];              %Volume fraction of grousers []
rocker_length = [];     %Length of rocker [m]
axle_width = [];        %Width of rocker [m]
rocker_height = [];     %Height between pivot of rocker and center of wheel when on flat terrain [m]

rover_fieldnames = {'fieldnames', 'r', 'b', 'h_g' 'r_s', 'W', 'mu', 'zeta', 'rocker_length', 'axle_width', 'rocker_height'};
rover = v2struct(rover_fieldnames);


%set soil properties (these values do not change thorughout terramechanics
% evaluations)
k = [];                 %Sinkage modulus [kN/m^(n+m+2)]
rho = [];               %Soil density [kg/m^3]
phi = [];               %Angle of repose [rad]
c = [];                 %Cohesion stress [kPa]
K = [];                 %Shear modulus [m]
X_c = [];               %Destructive angle of soil [rad]
n0 = [];                %Constant portion of slip dependent sinkage exponent []
n1 = [];                %Linear portion of slip dependent sinkage exponent []
n2 = [];                %Constant portion of skid dependent sinkage exponent []

soil_fieldnames = {'fieldnames', 'k', 'rho', 'phi', 'c', 'K', 'X_c', 'n0', 'n1', 'n2'};
soil = v2struct(soil_fieldnames);


%initialize terramechanics parameters (these values are changed throughout
%terramechanics evaluations)
theta_f = [];           %Entry angle of wheel [rad]
theta_r = [];           %Exit angle of wheel [rad]
theta_m = [];           %Angle of max normal stress under wheel [rad]
theta_0 = [];           %Angle of soil flow velocity sign transition under wheel [rad]
sigma = [];             %Normal stress under wheel at angle theta [kPa]
tau_l = [];             %Lateral shear stress under wheel at angle theta [kPa]
tau_t = [];             %Tangential shear stress under wheel at angle theta [kPa]
j_l = [];               %Lateral shear deformation module [m]
j_t = [];               %Tangential shear deformation module [m]
Rb = [];                %Bulldozing stress on side face of wheel [kPa]
h = [];                 %Sinkage of wheel [m]
sinkage_found = [];     %Whether the sinkage finding routine worked, if we used it
Fx = [];                %Computed value of force on wheel along wheel's x axis (forward)
Fy = [];                %Computed value of force on wheel along wheel's y axis (side face)
Fz = [];                %Computed value of force on wheel along wheel's z axis (up)


terr_fieldnames = {'fieldnames', 'theta_f', 'theta_r', 'theta_m', 'sigma', 'tau_l', 'tau_t', 'j_l', 'j_t', 'Rb', 'h', 'sinkage_found', 'Fx', 'Fy', 'Fz'};
terr = v2struct(terr_fieldnames);


%initialize tuned parameters (these values are not changed during
%individual terramechanics evaluations but may be optimized)
a0 = [];                %Empirical max stress angle constant offset []
a1 = [];                %Empirical max stress angle slip offset []
b0 = [];                %Empirical exit angle constant offset [deg]
b1 = [];                %Empirical exit angle slip offset [deg]
theta_m0 = [];          %Constant value of theta_m for a given slip angle for all skid values
theta_r0 = [];          %Constant value of theta_r for a given slip angle for all skid values
theta_f0 = [];          %Not actually tuned, this is the entrance angle for the given beta at slip=0

a_c = [];
a_b = [];
a_s = [];
a_sb = [];
b_c = [];
b_b = [];
b_s = [];
b_sb = [];

tuned_fieldnames = {'fieldnames', 'a0', 'a1', 'b0', 'b1', 'theta_m0', 'theta_r0', 'theta_f0', 'a_c', 'a_b', 'a_s', 'a_sb', 'b_c', 'b_b', 'b_s', 'b_sb'};
tuned = v2struct(tuned_fieldnames);


%initialize state parameters (these values do not change thorughout 
%terramechanics evaluations) but may be changed in a rover simulation
beta = [];              %Slip angle of wheel [rad]
slip = [];              %Slip ratio of wheel []
w = [];                 %Rotational (angular) velocity of wheel [rad/s]
v_x = [];               %Forward velocity of wheel [m/s]
v_y = [];               %Sideways velocity of wheel [m/s]
rover_w = 0;            %Yaw angular velocity of COM of rover [rad/s]
rover_v_x = [];          %Forward velocity of COM of rover [m/s]
rover_v_y = [];          %Sideways velocity of COM of rover [m/s]

state_fieldnames = {'fieldnames', 'beta', 'slip', 'w', 'v_x', 'v_y', 'rover_w', 'rover_v_x', 'rover_v_y'};
state = v2struct(state_fieldnames);

%initialize model options (these values should never change during model
%evaluation, they are initialized to default values so that they must be consciously set otherwise)
trench_on = 1;
ishigami_sidewall = 0;
linear_tuning_params = 0;

options_fieldnames = {'fieldnames', 'trench_on', 'ishigami_sidewall', 'linear_tuning_params'};
options = v2struct(options_fieldnames);


%make superstruct
fieldnames = {'fieldnames', 'rover', 'soil', 'terr', 'tuned', 'state', 'options'};
params = v2struct(fieldnames);
return