%This script is for tuning the values of a0, a1, b0, and b1 on real data
%collected in the terramechanics testbed. Run top subsections to choose the
%data you want to tune on. Save the tuned data as a .mat file to use later
%for modeling
clear all
warning('off', 'MATLAB:singularMatrix');%turn off singular matrix inversion warnings which we get from parallel lines in the profile intersection script
tuned_results = [];
Vrys = [0 3 5 8 10 12 20 40 100];
betas = [0 15 30 45 60 75 90];
% betas = [0]

exitflag = zeros(length(betas));

%Set starting values, data path, and rover parameters:
%Smooth wheel
%Starting guesses, x = [a0 a1 b0 b1 theta_m0 theta_r0 K c n0 n1 n2]
x0s = [0.37   0.50   -0.17   -0.01   0.26   -0.12    0.021    1.14    1.45    0.005    .75;... %0
       0.41   0.46   -0.18   -0.01   0.31   -0.14    0.021    1.14    1.45    0.005    .75;... %15
       0.60   0.21   -0.12   -0.15   0.52   -0.10    0.021    1.14    1.45    0.005    .75;... %30
       0.66   0.17   -0.14   -0.21   0.64   -0.14    0.021    1.14    1.45    0.005    .75;... %45
       0.76   0.12   -0.32   -0.10   0.84   -0.35    0.021    1.14    1.45    0.005    .75;... %60
       0.84   0.09   -0.48   -0.08   1.07   -0.61    0.021    1.14    1.45    0.005    .75;... %75
       0.91   0.07   -0.99   -0.01   1.53   -1.68    0.021    1.14    1.45    0.005    .75]'; %90
filepath = '../../terramech-testbed/nptm_control/process_data/all_smooth_data_2.mat'; %load smooth data
rovername = 'K10_mini_smooth'; %smooth wheel

% %Grousered wheel:
% %Starting guesses, x = [a0 a1 b0 b1 theta_m0 theta_r0 K c n0 n1 n2]
% x0s = [0.36   0.57   -0.72   -0.10    0.26   -0.51    0.0107    .10     1.45     0.05    .7;... %0
%        0.44   0.47   -0.59   -0.34    0.33   -0.45    0.021    1.14    1.45    0.005    .75;... %15
%        0.59   0.31   -0.45   -0.55    0.50   -0.39    0.021    1.14    1.45    0.005    .75;... %30
%        0.70   0.20   -0.44   -0.56    0.68   -0.43    0.021    1.14    1.45    0.005    .75;... %45
%        0.77   0.14   -0.44   -0.56    0.85   -0.49    0.021    1.14    1.45    0.005    .75;... %60
%        0.83   0.07   -0.40   -0.57    1.06   -0.51    0.021    1.14    1.45    0.005    .75;... %75
%        0.93   0.04   -0.40   -0.50    1.58   -0.51    0.021    1.14    1.45    0.005    .75]'; %90
% filepath = '../../terramech-testbed/nptm_control/process_data/all_grouser_data_2.mat'; %load grouser data
% rovername = 'K10_mini'; %Grousered wheel

%Tune over both slip and skid conditions for each slip angle separately
for i=1:length(betas)
    % Load Wheel Data
    clearvars -except i betas Vrys filepath tuned_results x0s rovername
    load(filepath); %load grouser data
    x0 = x0s(:,i);

    %Average the data and extract just the slip angle we want
    [Fx_data, Fy_data, Fz_data, Z_data, slip_data, V_data, Vry_data, angle_data, ~, ~, ~, ~, ~,~] = average_data_to_model(all_results, betas, Vrys);
    Fx_data = Fx_data(:,i);
    Fy_data = Fy_data(:,i);
    Fz_data = Fz_data(:,i);
    Z_data = Z_data(:,i);
    h0 = Z_data(5)/1000; %sinkage at slip=0, changing Vrys vector will break this
    min_h_skid = min(Z_data(1:5))/1000;
    slip_data = slip_data(:,i);
    V_data = V_data(:,i);
    Vry_data = Vry_data(:,i);
    angle_data = angle_data(:,i);
    Vry = Vrys;
    beta = betas(i);
    
    %repackage the data into a struct so the optimization can just use it
    %without modification
    all_results = [];
    for j=1:length(Vrys)
        result.avg_Fx = -Fy_data(j);
        result.avg_Fy = -Fx_data(j);
        result.avg_Fz = -Fz_data(j);
        result.avg_Z = -Z_data(j);
        result.slip = slip_data(j);
        result.avg_Vx = V_data(j);
        result.avg_Vry = Vry_data(j);
        result.avg_angle = angle_data(j);
        result.beta = beta;
        result.Vry = Vrys(j);
        all_results = [all_results, result];
    end    
    
    [x, fval, exitflag(i)] = a_b_tuning_optimization(all_results, h0, min_h_skid, x0, rovername);
    fprintf('beta = %d, exitflag = %d\n', betas(i), exitflag(i))
    
    % Save the tuned values of a and b
    for k=1:length(all_results)
        all_results(k).a0 = x(1);
        all_results(k).a1 = x(2);
        all_results(k).b0 = x(3);
        all_results(k).b1 = x(4);
        all_results(k).theta_m0 = [];
        all_results(k).theta_r0 = [];
        all_results(k).K = x(7);
        all_results(k).c = x(8);
        all_results(k).n0 = x(9);
        all_results(k).n1 = x(10);
        all_results(k).n2 = x(11);
        
        
    end
    tuned_results = [tuned_results; all_results'];
end
warning('on', 'MATLAB:singularMatrix'); %put the warnings back on
