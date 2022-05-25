% Set the rover geometry values needed to run terramechanics simulations. 
% User must specify rover name from the following list:  K10_mini, 
% K10_mini_smooth, WD2, LATUV, or KRex

function vars = set_rover_params(params, rovername)
if strcmp(rovername,'K10_mini_grouser') %This is actually values for the adjustable rover with K10 mini's wheels
    params.rover.r = .058; %Wheel radius
    params.rover.b = .05; %Wheel width
    params.rover.W = 10/1000; %Based on adjustable rover weight of 4kg, so 1kg on each wheel = 10N on each wheel
    params.rover.mu = .1; %Area ratio of grouser tips
    params.rover.rocker_length = .3; %Rocker length (between wheel centers) in m
    params.rover.axle_width = .3; %Axle width (between wheel centers) in m
    params.rover.rocker_height = .09202; %height between pivot of rocker and center of wheel when on flat terrain
    params.rover.h_g = 5/1000; %grouser height in m
    params.rover.r_s = params.rover.r + params.rover.h_g; %shearing radius for soil under wheel in m, guesstimated
    params.rover.zeta = .2; %volume fraction of grousers for trench soil transport
    
elseif strcmp(rovername,'K10_mini_smooth') %This is actually values for the adjustable rover with K10 mini's wheels, without grousers
    params.rover.r = .058; %Wheel radius
    params.rover.b = .05; %Wheel width
    params.rover.W = 10/1000; %Based on adjustable rover weight of 4kg, so 1kg on each wheel = 10N on each wheel
    params.rover.mu = 1; %Area ratio of grouser tips
    params.rover.rocker_length = .3; %Rocker length (between wheel centers) in m
    params.rover.axle_width = .3; %Axle width (between wheel centers) in m
    params.rover.rocker_height = .09202; %height between pivot of rocker and center of wheel when on flat terrain
    params.rover.h_g = 0; %grouser height in m
    params.rover.r_s = params.rover.r + params.rover.h_g; %shearing radius for soil under wheel in m, guesstimated
    params.rover.r_s = params.rover.r*1.02; 
    params.rover.zeta = .1; %volume fraction of grousers for trench soil transport
    fprintf('r_s = 1.02r\n')
    
elseif strcmp(rovername,'WD2')
    params.rover.r = .25; %[m]
    params.rover.b = .32; %[m]
    params.rover.W = (65.7)*10/1000; %[?]
    params.rover.mu = .1; %guesstimated
    params.rover.rocker_length = [];
    params.rover.axle_width = [];
    params.rover.rocker_height = [];
    params.rover.h_g = .075/100;
    params.rover.rs = .075/100+params.rover.r; %[m]
    params.rover.zeta = 1; %volume fraction of grousers for trench soil transport
    %params.state.vx = 3.2/100; %[m/s]
    %params.state.w = (4/100)/params.rover.r; %[rad/sec]
    %params.state.s = 0.2;
    %params.state.beta = 0;
    
    
    %{
elseif strcmp(rovername,'MER')
    %estimated values
    print('MER values not available yet. Please rerun with {default, K10_mini, LATUV, KRex} as rover_name')
    params.r = .26/2;
    params.b = .1;
    params.W = ;
    params.mu = ;
    params.rocker_length = ;
    params.axle_width = ;
    params.pivot_height = ;
    params.rs = ;
    %}
    
elseif strcmp(rovername, 'LATUV')
    %These values are estimated, should get real #s from Dimi
    params.rover.r = .2;
    params.rover.b = .2;
    params.rover.W = 3060/1000/4; %Weight in kN from STTR final report
    params.rover.mu = .5;
    params.rover.h_g = .1*params.rover.r; 
    params.rover.rocker_length = 1;
    params.rover.axle_width = 1;
    params.rover.rocker_height = .25;
    params.rover.r_s = 1.1*params.rover.r;
    params.rover.zeta = 1; %volume fraction of grousers for trench soil transport
    
elseif strcmp(rovername, 'KRex')
    params.rover.r = .46/2; %Measured to area of ground contact in the Atacama
    params.rover.b = .152; %from STTR final presentation
    params.rover.W = 2290/1000/4; %Weight in kN from STTR final report
    params.rover.mu = .1; %guesstimated
    params.rover.rocker_length = 1.4;%from STTR final presentation
    params.rover.axle_width = 1.15 ; %from STTR final presentation
    params.rover.rocker_height = (.64+.42)/2 - .508/2; %from STTR final presentation
    params.rover.r_s = 1.1*params.rover.r; %guesstimated
    params.rover.h_g = 0.015; %guesstimated
    params.rover.zeta = 1; %volume fraction of grousers for trench soil transport
    params.terr.phi = deg2rad(33.5); %average value measured by inclinometer in the atacama desert
    
    
% elseif strcmp(rovername, 'smooth')
%     params.rover.r = .048; %Wheel radius
%     params.rover.b = .05; %Wheel width
%     params.rover.W = 10/1000; %I don't know the actual rover weight
%     params.rover.mu = 0; %Area ratio of grouser tips
%     params.rover.rocker_length = .3; %Rocker length (between wheel centers) in m
%     params.rover.axle_width = .3; %Axle width (between wheel centers) in m
%     params.rover.rocker_height = .09202; %height between pivot of rocker and center of wheel when on flat terrain
%     params.rover.r_s = 1.1*params.r; %shearing radius for soil under wheel in m, guesstimated
%     params.rover.zeta = 1; %volume fraction of grousers for trench soil transport
%     params.rover.h_g = 0.0;
    
% elseif strcmp(rovername, 'JiaWheel')
%     params.rover.r = 0.092;
%     params.rover.b = 0.107;
%     params.rover.W = 64/1000;
%     params.rover.mu = 1;
%     params.rover.r_s = 0.092;
%     params.rover.h_g = 0.001;
%     params.rover.zeta = 1; %volume fraction of grousers for trench soil transport
else
    error('Rover name not recognized. Valid rover names: K10_mini_grouser, K10_mini_smooth, WD2, LATUV, or KRex.\n')
end
vars = params;
end