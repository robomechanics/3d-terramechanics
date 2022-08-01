%% Load already processed data for the smooth wheel
load('all_grouser_data.mat'); %load grouser data
%% Load already processed data for the smooth wheel
load('all_smooth_data.mat'); %load smooth data
%% Process the raw data by extracting average values from the steady state portion of each trial
min_time = 70;
max_time = 75;
disp_plots = 0; %Change to 1 to display raw data for each trial as a time series plot
def_fit = testbed_deflection_analysis(0); %Get linear compensation for vertical deflection of the testbed under load
def_fit = fit([1 0]',[0 0]', 'poly1');

% foldernames = {'gr810_vx_10/11-16-20/'; ... %Data for smooth wheel
%                'gr630_vx_10/11-19-20/'; ...
%                'gr630_vx_10/11-23-20/'; ...
%                'gr810_vx_10/11-24-20/'; ...
%                'gr810_vx_10/11-25-20/'; ...
%                'gr810_vx_10/11-26-20/'; ...
%                'gr810_vx_10/11-27-20/'; ...
%                'gr810_vx_10/11-30-20/'; ...
%                'gr810_vx_10/12-2-20/'; ...
%                'gr630_vx_10/12-3-20/'; ...
%                'gr810_vx_10/12-5-20/'; ...
%                'gr810_vx_10/12-7-20/'; ...
%                'gr810_vx_10/12-8-20/'; ...
%                'gr810_vx_10/12-10-20/'; ...
%                };
% DST = zeros(14,1);
           
foldernames = {'gr810_grouser/12-14-20/'; ... %Data for grousered wheel           
               'gr810_grouser/3-2-21/'; ...
               'gr810_grouser/3-3-21/'; ...
               'gr810_grouser/3-4-21/'; ...
               'gr810_grouser/3-9-21/'; ...
               'gr630_grouser/3-11-21/'; ...
               'gr630_grouser/3-16-21/'; ...
    }; 
DST = [0 0 0 0 0 0 1];

%Each entry is in the form of [Vry, beta aka angle, trialnum]. Note that
%trialnum is the last digit in the filename, and will usually be 1 or 2 as
%it reflects how many trials for that set of (Vry, beta) occurred on that
%day.
% test_speed_angle_n = { %Data for smooth wheel
%                         {[3,30,1]; [3,60,1]; [5,0,1]; [8,0,1]; [10,30,1]; [10,60,1]; [12,0,1]; [12,30,1]; [20,60,1]}; ... %11-16-20 (0-10)
%                         {[40,0,1]; [40,15,1]; [40,30,1]; [40,45,1]; [40,60,1]; [40,75,1]; [40,90,1]; [100,75,1]}; ... %11-19-20 (11-18)
%                         {[40,0,1]; [40,15,1]; [40,30,1]; [40,45,1]; [40,60,1]; [40,75,1]; [40,90,1]; [100,15,1]; [100,30,1]; [100,45,1]; [100,60,1]}; ... %11-23-20 (19-29)
%                         {[0,15,1]}; ... %11-24-20 (30)
%                         {[0,0,1]}; ... %11-25-20 (31)
%                         {[0,30,1]; [0,45,1]; [0,75,1]; [0,90,1]}; %11-26-20 (32-35)
%                         {[0,15,1]; [0,45,1]; [0,60,1]; [0,90,1]}; %11-27-20 (36-39)
%                         {[0,0,1]; [0,30,1]; [0,60,1]; [0,75,1]}; %11-30-20 (40-43)
%                         {[0,15,1]; [0,45,1]; [0,60,1]; [0,90,1]}; %12-2-20 (44-47)
%                         {[100,0,1]; [100,0,2]; [100,0,3]; [100,15,1]; [100,15,2]; [100,30,1]; [100,30,2]; [100,45,1]; ... %12-3-20 (48-55)
%                            [100,45,2]; [100,60,1]; [100,60,2]; [100,75,1]; [100,75,2]; [100,90,1]; [100,90,2]; [100,90,3]}; %12-3-20 (56-63)
%                         {[3,0,1]; [3,15,1]; [3,30,1]; [3,45,1]; [3,75,1]; [3,90,1]; [5,15,1]; [5,30,1]; [5,45,1]; [5,60,1]; [5,75,1]; [5,90,1]; ... %12-5-20 (64-75)
%                            [8,15,1]; [8,15,2]; [8,30,1]; [8,45,1]; [8,60,1]; [8,75,1]; [8,90,1]; ... %12-5-20 (76-82)
%                            [10,0,1]; [10,15,1]; [10,30,1]; [10,45,1]; [10,75,1]; [10,75,2]; [10,90,1]; ... %12-5-20 (83-89)
%                            [20,0,1]; [20,15,1]; [20,30,1]; [20,45,1]; [20,75,1]; [20,90,1]; [20,90,2]}; ... %12-5-20 (90-96)
%                         {[3,0,1]; [3,0,2]; [3,15,1]; [3,30,1]; [3,45,1]; [3,60,1]; [3,75,1]; [3,75,2]; [3,90,1]; [3,90,2]; ... %12-7-20 (97-106)
%                             [5,0,1]; [5,15,1]; [5,15,2]; [5,30,1]; [5,30,2]; [5,45,1]; [5,60,1]; ... %12-7-20 (107-113)
%                             [8,0,1]; [8,0,2]; [8,15,1]; [8,30,1]; [8,45,1]; [8,60,1]; [8,60,2]; [8,75,1]; [8,90,1]; [8,90,2]; ... %12-7-20 (114-124)
%                             [10,0,1]; [10,0,2]; [10,15,1]; [10,30,1]; [10,45,1]; [10,60,1]; [10,60,2]; [10,75,1]; [10,90,1]; [10,90,2]; ... %12-7-20 (125-134)
%                             [12,45,1]; [20,0,1]; [20,15,1]; [20,15,2]; [20,30,1]; [20,45,1]; [20,45,2]; [20,60,1]; [20,60,2]; [20,75,1]; [20,75,2]; ... %12-7-20 (135-145)
%                             [40,0,1]; [40,30,1]; [40,45,1]; [40,90,1]}; ... %12-7-20 (146-149)
%                         {[0,0,1]; [0,30,1]; [0,75,1]; [3,15,1]; [3,45,1]; [3,60,1]; [5,0,1]; [5,75,1]; [5,90,1]; [8,30,1]; [8,45,1]; [8,75,1]; ... %12-8-20 (150-160)
%                             [10,15,1]; [10,45,1]; [12,30,1]; [12,60,1]; [12,90,1]; [20,0,1]; [20,30,1]; [20,90,1]; [40,15,1]; [40,60,1]; [40,75,1]}; ... %12-8-20 (161-171)
%                         {[5,45,1]; [5,60,1]; [5,75,1]; [5,90,1]; [12,0,1]; [12,0,2]; [12,15,1]; [12,15,2]; [12,15,3]; [12,30,1]; ... (172-181) %12-10-20
%                            [12,45,1]; [12,45,2]; [12,60,1]; [12,60,2]; [12,75,1]; [12,75,2]; [12,75,3]; [12,90,1]; [12,90,2]} %12-10-20 (182-190)
%                      };

test_speed_angle_n = { %Data for grousered wheel
                        {[0,0,1]; [5,0,1]; [10,0,1]; [20,0,1]; [40,0,1]; [0,30,1]; [5,30,1]; [10,30,1]; [20,30,1]; [40,30,1]; ... %12-14-20
                        [0,60,1]; [5,60,1]; [10,60,1]; [20,60,1]; [40,60,1]; [0,90,1]; [5,90,1]; [10,90,1]; [20,90,1]; [40,90,1]}; ... %12-14-20
                        {[0,15,1]; [0,45,1]; [0,75,1]; [0,90,1]; [3,0,1]; [3,15,1]; [3,30,1]; [3,45,1]; [3,60,1]; [3,75,1]; [3,90,1]; ... %3-2-21
                        [5,15,1]; [5,45,1]; [5,60,1]; [5,75,1]; [8,0,1]; [8,15,1]; [8,30,1]; [8,45,1]; [8,60,1]; [8,75,1]; [8,90,1]; ... %3-2-21
                        [10,15,1]; [10,30,1]; [10,45,1]; [10,75,1]; [12,15,1]; [12,30,1]; [12,45,1]; [12,60,1]; [12,75,1]; [12,90,1]; ... %3-2-21
                        [20,0,1]; [20,15,1]; [20,45,1]; [20,75,1]; [40,15,1]; [40,45,1]; [40,75,1]; [12,0,1]}; ... % 3-2-21
                        {[0,0,1]; [0,30,1]; [0,60,1]; [0,75,1]; [3,0,1]; [3,15,1]; [3,45,1]; [3,75,1]; [3,90,1]; [5,15,1]; [5,30,1]; ... %3-3-21
                        [5,45,1]; [5,90,1]; [8,0,1]; [8,15,1]; [8,30,1]; [8,45,1]; [8,60,1]; [8,75,1]; [10,0,1]; [10,15,1]; [10,45,1]; ... %3-3-21
                        [10,60,1]; [10,90,1]; [12,0,1]; [12,15,1]; [12,30,1]; [12,60,1]; [12,75,1]; [20,30,1]; [20,45,1]; [20,60,1]; [20,75,1]; [20,90,1]}; ... %3-3-21
                        {[0,15,1]; [0,45,1]; [0,45,2]; [0,60,1]; [3,30,1]; [3,30,2]; [3,60,1]; [3,90,1]; [5,0,1]; [5,15,1]; [5,30,1]; [5,75,1]; ... %3-4-21
                        [8,45,1]; [8,90,1]; [8,90,2]; [10,0,1]; [10,75,1]; [10,75,2]; [12,15,1]; [12,45,1]; [12,75,1]; [12,90,1]; [20,0,1]; [20,15,1]; [20,60,1]}; ... %3-4-21
                        {[0,0,1]; [0,15,1]; [0,30,1]; [0,75,1]; [0,90,1]; [3,0,1]; [3,15,1]; [3,45,1]; [3,60,1]; [3,75,1]; ... %3-9-21
                        [5,0,1]; [5,45,1]; [5,60,1]; [5,75,1]; [5,90,1]; [8,0,1]; [8,15,1]; [8,30,1]; [8,60,1]; [8,75,1]; [10, 15,1]; [10,30,1]; [10,45,1]; [10,60,1]; ... %3-9-21
                        [10,90,1]; [12,0,1]; [12,30,1]; [12,45,1]; [12,60,1]; [12,90,1]; [20,15,1]; [20,30,1]; [20,45,1]; [20,75,1]; [20,90,1]}; ... %3-9-21
                        {[40,15,1]; [40,30,1]; [40,60,1]; [40,90,1]; [100,0,1]; [100,30,1]; [100,45,1]; [100,75,1]}; ... %3-11-21 
                        {[40,0,1]; [40,0,2]; [40,15,1]; [40,30,1]; [40,45,1]; [40,45,2]; ... %3-16-21
                        [40,60,1]; [40,75,1]; [40,75,2]; [40,90,1]; [100,0,1]; [100,0,2]; [100,15,1]; [100,15,2]; [100,15,3]; ... %3-16-21
                        [100,30,1]; [100,30,2]; [100,45,1]; [100,45,2]; [100,60,1]; [100,60,2]; [100,60,3]; [100,75,1]; [100,75,2]; [100,90,1]; [100,90,2]; [100,90,3]} ... %3-16-21
                    };

all_results = [];

for i=1:length(foldernames)
    files = test_speed_angle_n{i};
    foldername = foldernames{i};
    for j=1:length(files)
        val = files{j};
        filename = sprintf('10_%.0f_%.0f_%i', val(1), val(2), val(3));
        results = read_ft_speed_data(foldername, filename, min_time, max_time, disp_plots, DST(i), def_fit);
        results.Vry = val(1);
        results.beta = val(2);
        results.trial_num = val(3);
        if val(1) == 0
            results.avg_Vry = 0;
            results.std_Vry = 0;
        end
                       
        results.slip = (results.avg_Vry - 10)/max(10, results.avg_Vry);
        all_results = [all_results;results];      
        
    end
end
%% Forces Figure
figure()
sgtitle ('Wheel Forces')
subplot(1,3,2)
title('Fy (Sidewall)')
axis([-1.1 1 -25 50])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,3,1)
title('Fx (Tractive)')
axis([-1.1 1 -25 50])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,3,3)
title('Fz (Load)')
axis([-1.1 1 -25 50])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

for i=1:length(all_results)
    result = all_results(i);
    %Select color of plotted point
    switch result.beta
        case 0
            color = cmuColor('red-web');
        case 15
            color = cmuColor('gold');
        case 30
            color = cmuColor('teal');
        case 45
            color = cmuColor('sky-blue');
        case 60
            color = cmuColor('palladian-green');
        case 75
            color = cmuColor('blue-thread');
        case 90
            color = cmuColor('scots-rose');
        otherwise
            color = 'k';
    end
    
    subplot(1,3,2)
    plot(result.slip, -result.avg_Fx, 'o', 'MarkerEdgeColor', color)

    subplot(1,3,1)
    plot(result.slip, -result.avg_Fy, 'o', 'MarkerEdgeColor', color)
    
    subplot(1,3,3)
    plot(result.slip, -result.avg_Fz, 'o', 'MarkerEdgeColor', color)
    
end
leg(1) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('red-web'));
leg(2) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('gold'));
leg(3) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('teal'));
leg(4) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('sky-blue'));
leg(5) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('palladian-green'));
leg(6) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('blue-thread'));
leg(7) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('scots-rose')); 
legg = legend(leg, '0', '15', '30', '45', '60', '75', '90');
title(legg, 'Slip Angle');
hold off

%% Sinkage Figure
figure()
sgtitle('Sinkage')
xlabel('Slip')
ylabel('Depth (mm)')
hold on
for i=1:length(all_results)
    result = all_results(i);
    %Select color of plotted point
    switch result.beta
        case 0
            color = cmuColor('red-web');
        case 15
            color = cmuColor('gold');
        case 30
            color = cmuColor('teal');
        case 45
            color = cmuColor('sky-blue');
        case 60
            color = cmuColor('palladian-green');
        case 75
            color = cmuColor('blue-thread');
        case 90
            color = cmuColor('scots-rose');
        otherwise
            color = 'k';
    end
    
    plot(result.slip, -result.avg_Z, 'o', 'MarkerEdgeColor', color)
    
end
leg(1) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('red-web'));
leg(2) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('gold'));
leg(3) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('teal'));
leg(4) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('sky-blue'));
leg(5) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('palladian-green'));
leg(6) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('blue-thread'));
leg(7) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('scots-rose')); 
legg = legend(leg, '0', '15', '30', '45', '60', '75', '90', 'Location', 'northeastoutside');
title(legg, 'Slip Angle');
hold off

%% Slope coefficients Figure
figure()
sgtitle ('Wheel Forces')
subplot(2,3,2)
title('Fy (Sidewall)')
axis([-1.1 1 -.5 1])
hold on
xlabel('Slip Ratio')
ylabel('Slope (N/s)')

subplot(2,3,1)
title('Fx (Tractive)')
axis([-1.1 1 -.5 1])
hold on
xlabel('Slip Ratio')
ylabel('Slope (N/s)')

subplot(2,3,3)
title('Fz (Load)')
axis([-1.1 1 -.5 1])
hold on
xlabel('Slip Ratio')
ylabel('Slope (N/s)')

subplot(2,3,5)
title('Fy (Sidewall)')
axis([-1.1 1 -30 50])
hold on
xlabel('Slip Ratio')
ylabel('Offset (N)')

subplot(2,3,4)
title('Fx (Tractive)')
axis([-1.1 1 -30 50])
hold on
xlabel('Slip Ratio')
ylabel('Offset (N)')

subplot(2,3,6)
title('Fz (Load)')
axis([-1.1 1 -30 50])
hold on
xlabel('Slip Ratio')
ylabel('Offset (N)')

for i=1:length(all_results)
    result = all_results(i);
    %Select color of plotted point
    switch result.beta
        case 0
            color = cmuColor('red-web');
        case 15
            color = cmuColor('gold');
        case 30
            color = cmuColor('teal');
        case 45
            color = cmuColor('sky-blue');
        case 60
            color = cmuColor('palladian-green');
        case 75
            color = cmuColor('blue-thread');
        case 90
            color = cmuColor('scots-rose');
        otherwise
            color = 'k';
    end
    
    subplot(2,3,2)
    plot(result.slip, result.polyx(1), 'o', 'MarkerEdgeColor', color)

    subplot(2,3,1)
    plot(result.slip, -result.polyy(1), 'o', 'MarkerEdgeColor', color)
    
    subplot(2,3,3)
    plot(result.slip, -result.polyz(1), 'o', 'MarkerEdgeColor', color)
    
    subplot(2,3,5)
    plot(result.slip, result.polyx(2), 'o', 'MarkerEdgeColor', color)

    subplot(2,3,4)
    plot(result.slip, -result.polyy(2), 'o', 'MarkerEdgeColor', color)
    
    subplot(2,3,6)
    plot(result.slip, -result.polyz(2), 'o', 'MarkerEdgeColor', color)
    
end
leg(1) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('red-web'));
leg(2) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('gold'));
leg(3) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('teal'));
leg(4) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('sky-blue'));
leg(5) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('palladian-green'));
leg(6) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('blue-thread'));
leg(7) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('scots-rose')); 
legend(leg, '0', '15', '30', '45', '60', '75', '90')
hold off

%% Data variance values/plots, highly inefficient

betas = [0 15 30 45 60 75 90];
Vrys = [0 3 5 8 10 12 20 40 100];

Fx_means = zeros(length(betas), length(Vrys));
Fx_stds = zeros(length(betas), length(Vrys));
Fy_means = zeros(length(betas), length(Vrys));
Fy_stds = zeros(length(betas), length(Vrys));
Fz_means = zeros(length(betas), length(Vrys));
Fz_stds = zeros(length(betas), length(Vrys));
Z_means = zeros(length(betas), length(Vrys));
Z_stds = zeros(length(betas), length(Vrys));
slip_means = zeros(length(betas), length(Vrys));
slip_stds = zeros(length(betas), length(Vrys));

temp = [];
Fx = [];
Fy = [];
Fz = [];
Z = [];
slip = [];

for j=1:length(betas)
    for k=1:length(Vrys)
        for i=1:length(all_results)
            datapt = all_results(i);
            if (datapt.beta == betas(j) && datapt.Vry == Vrys(k))
                temp = [temp, datapt];
            end
        end
        for m=1:length(temp)
            datapt = temp(m);
            Fx = [Fx, -datapt.avg_Fy];
            Fy = [Fy, -datapt.avg_Fx];
            Fz = [Fz, -datapt.avg_Fz];
            Z = [Z, -datapt.avg_Z];
            slip = [slip, datapt.slip];
            
        end
        Fx_means(j,k) = mean(Fx);
        Fx_stds(j,k) = std(Fx);
        Fy_means(j,k) = mean(Fy);
        Fy_stds(j,k) = std(Fy);
        Fz_means(j,k) = mean(Fz);
        Fz_stds(j,k) = std(Fz);
        Z_means(j,k) = mean(Z);
        Z_stds(j,k) = std(Z);
        slip_means(j,k) = mean(slip);
        slip_stds(j,k) = std(slip);
        
        temp = [];
        Fx = [];
        Fy = [];
        Fz = [];
        Z = [];
        slip = [];
    end
end
    
Fx_variation = Fx_stds./Fx_means;
Fy_variation = Fy_stds./Fy_means;
Fz_variation = Fz_stds./Fz_means;
Z_variation = Z_stds./Z_means;

