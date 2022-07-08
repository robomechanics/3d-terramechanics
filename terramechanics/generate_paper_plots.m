% Generates all plots used in JTM paper by running the model with various
% options.

%Set up for all versions
clear all
warning('off', 'MATLAB:singularMatrix')
betas = [0 15 30 45 60 75 90];
Vrys = [0 3 5 8 10 12 20 40 100];

%Set up options for various versions
wheel_g = 'K10_mini_grouser';
wheel_s = 'K10_mini_smooth';

%Load and average grousered data
load('trench_off_grouser.mat')
all_results_g = tuned_results;
[Fx_data_g, Fy_data_g, Fz_data_g, Z_data_g, slip_data_g, V_data_g, Vry_data_g, angle_data_g, a0s_g, a1s_g, b0s_g, b1s_g, theta_m0s_g, theta_r0s_g] = average_data_to_model(all_results_g, betas, Vrys);

%Load and average smooth data
load('trench_off_smooth.mat')
all_results_s = tuned_results;
[Fx_data_s, Fy_data_s, Fz_data_s, Z_data_s, slip_data_s, V_data_s, Vry_data_s, angle_data_s, a0s_s, a1s_s, b0s_s, b1s_s, theta_m0s_s, theta_r0s_s] = average_data_to_model(all_results_s, betas, Vrys);

%Run model with grousers 
[all_params_gi, Fx_approx_gi, Fy_approx_gi, h_gi, sinkage_found_gi] = run_model(wheel_g, betas, Vrys, 1, 1, Fz_data_g, V_data_g, Vry_data_g, a0s_g, a1s_g, b0s_g, b1s_g); %Grousered, Ishigami
[all_params_gn, Fx_approx_gn, Fy_approx_gn, h_gn, sinkage_found_gn] = run_model(wheel_g, betas, Vrys, 0, 0, Fz_data_g, V_data_g, Vry_data_g, a0s_g, a1s_g, b0s_g, b1s_g); %Grousered, no trench
[all_params_gt, Fx_approx_gt, Fy_approx_gt, h_gt, sinkage_found_gt] = run_model(wheel_g, betas, Vrys, 1, 0, Fz_data_g, V_data_g, Vry_data_g, a0s_g, a1s_g, b0s_g, b1s_g); %Grousered, trench

%Run model with smooth wheel 
[all_params_si, Fx_approx_si, Fy_approx_si, h_si, sinkage_found_si] = run_model(wheel_s, betas, Vrys, 1, 1, Fz_data_s, V_data_s, Vry_data_s, a0s_s, a1s_s, b0s_s, b1s_s); %Smooth, Ishigami
[all_params_sn, Fx_approx_sn, Fy_approx_sn, h_sn, sinkage_found_sn] = run_model(wheel_s, betas, Vrys, 0, 0, Fz_data_s, V_data_s, Vry_data_s, a0s_s, a1s_s, b0s_s, b1s_s); %Smooth, no trench
[all_params_st, Fx_approx_st, Fy_approx_st, h_st, sinkage_found_st] = run_model(wheel_s, betas, Vrys, 1, 0, Fz_data_s, V_data_s, Vry_data_s, a0s_s, a1s_s, b0s_s, b1s_s); %Smooth, trench

%Calculate errors on grouser model
[avg_Fx_error_norm_gn, std_Fx_error_norm_gn, avg_Fy_error_norm_gn, std_Fy_error_norm_gn, avg_Z_error_norm_gn, std_Z_error_norm_gn] = model_errors(Fx_data_g, Fy_data_g, Fz_data_g, Z_data_g, Fx_approx_gn, Fy_approx_gn, h_gn, all_params_gn(1));
[avg_Fx_error_norm_gt, std_Fx_error_norm_gt, avg_Fy_error_norm_gt, std_Fy_error_norm_gt, avg_Z_error_norm_gt, std_Z_error_norm_gt] = model_errors(Fx_data_g, Fy_data_g, Fz_data_g, Z_data_g, Fx_approx_gt, Fy_approx_gt, h_gt, all_params_gt(1));

%Calculate errors on smooth model
[avg_Fx_error_norm_sn, std_Fx_error_norm_sn, avg_Fy_error_norm_sn, std_Fy_error_norm_sn, avg_Z_error_norm_sn, std_Z_error_norm_sn] = model_errors(Fx_data_s, Fy_data_s, Fz_data_s, Z_data_s, Fx_approx_sn, Fy_approx_sn, h_sn, all_params_sn(1));
[avg_Fx_error_norm_st, std_Fx_error_norm_st, avg_Fy_error_norm_st, std_Fy_error_norm_st, avg_Z_error_norm_st, std_Z_error_norm_st] = model_errors(Fx_data_s, Fy_data_s, Fz_data_s, Z_data_s, Fx_approx_st, Fy_approx_st, h_st, all_params_st(1));

%Generate plots and tables
fprintf('Grouser, no trench:\n')
print_tables(all_params_gn, betas, avg_Fx_error_norm_gn, std_Fx_error_norm_gn, avg_Fy_error_norm_gn, std_Fy_error_norm_gn, avg_Z_error_norm_gn, std_Z_error_norm_gn, a0s_g, a1s_g, b0s_g, b1s_g); %Grouser, no trench
fprintf('Grouser, trench:\n')
print_tables(all_params_gt, betas, avg_Fx_error_norm_gt, std_Fx_error_norm_gt, avg_Fy_error_norm_gt, std_Fy_error_norm_gt, avg_Z_error_norm_gt, std_Z_error_norm_gt, a0s_g, a1s_g, b0s_g, b1s_g); %Grouser, trench
fprintf('Smooth, no trench:\n')
print_tables(all_params_sn, betas, avg_Fx_error_norm_sn, std_Fx_error_norm_sn, avg_Fy_error_norm_sn, std_Fy_error_norm_sn, avg_Z_error_norm_sn, std_Z_error_norm_sn, a0s_s, a1s_s, b0s_s, b1s_s); %Smooth, no trench
fprintf('Smooth, trench:\n')
print_tables(all_params_st, betas, avg_Fx_error_norm_st, std_Fx_error_norm_st, avg_Fy_error_norm_st, std_Fy_error_norm_st, avg_Z_error_norm_st, std_Z_error_norm_st, a0s_s, a1s_s, b0s_s, b1s_s); %Smooth, trench

%Grousered force and sinkage plot
[fig_g, ax_g, tiles_g] = force_sinkage_plot(all_params_gi, betas, Vrys, slip_data_g, Fx_data_g, Fy_data_g, Z_data_g, Fx_approx_gi, Fy_approx_gi, h_gi, sinkage_found_gi); %Grousered, Ishigami
force_sinkage_plot(all_params_gn, betas, Vrys, slip_data_g, Fx_data_g, Fy_data_g, Z_data_g, Fx_approx_gn, Fy_approx_gn, h_gn, sinkage_found_gn, fig_g, ax_g, tiles_g); %Grousered, no trench
force_sinkage_plot(all_params_gt, betas, Vrys, slip_data_g, Fx_data_g, Fy_data_g, Z_data_g, Fx_approx_gt, Fy_approx_gt, h_gt, sinkage_found_gt, fig_g, ax_g, tiles_g); %Grousered, trench

%Smooth force and sinkage plot
[fig_s, ax_s, tiles_s] = force_sinkage_plot(all_params_si, betas, Vrys, slip_data_s, Fx_data_s, Fy_data_s, Z_data_s, Fx_approx_si, Fy_approx_si, h_si, sinkage_found_si); %Smooth, Ishigami
force_sinkage_plot(all_params_sn, betas, Vrys, slip_data_s, Fx_data_s, Fy_data_s, Z_data_s, Fx_approx_sn, Fy_approx_sn, h_sn, sinkage_found_sn, fig_s, ax_s, tiles_s); %Smooth, no trench
force_sinkage_plot(all_params_st, betas, Vrys, slip_data_s, Fx_data_s, Fy_data_s, Z_data_s, Fx_approx_st, Fy_approx_st, h_st, sinkage_found_st, fig_s, ax_s, tiles_s); %Smooth, trench

%Legend plot
legend_plot();

%Grousered tuning parameters plot
tuning_params_plot(all_params_gn, betas, a0s_g, a1s_g); %Grousered, no trench

%Smooth tuning parameters plot
tuning_params_plot(all_params_sn, betas, a0s_s, a1s_s); %Smooth, no trench

warning('on', 'MATLAB:singularMatrix')

function [avg_Fx_error_norm, std_Fx_error_norm, avg_Fy_error_norm, std_Fy_error_norm, avg_Z_error_norm, std_Z_error_norm] = model_errors(Fx_data, Fy_data, Fz_data, Z_data, Fx_approx, Fy_approx, h, params)
%Calculates error statistics for the model generated forces and sinkages
r = params.rover.r;
Fx_error = 1000*Fx_approx(1:end,:) - Fx_data(1:end,:);
Fy_error = 1000*Fy_approx(1:end,:) - Fy_data(1:end,:);
%Fz_error = 1000*Fz_approx(1:end,:) - Fz_data(1:end,:);
Z_error = 1000*h(1:end,:) - Z_data(1:end,:);
avg_Fx_error_norm = mean(Fx_error./Fz_data(1:end,:)*100, 1);
%med_Fx_error_norm = median(Fx_error./Fz_data(1:end,:)*100, 1);
std_Fx_error_norm = std(Fx_error./Fz_data(1:end,:)*100, 1);
avg_Fy_error_norm = mean(Fy_error./Fz_data(1:end,:)*100, 1);
%med_Fy_error_norm = median(Fy_error./Fz_data(1:end,:)*100, 1);
std_Fy_error_norm = std(Fy_error./Fz_data(1:end,:)*100, 1);
avg_Z_error_norm = mean(Z_error./r/1000*100, 1);
%med_Z_error_norm = median(Z_error./r/1000*100, 1);
std_Z_error_norm = std(Z_error./r/1000*100, 1);
end

function print_tables(all_params, betas, avg_Fx_error_norm, std_Fx_error_norm, avg_Fy_error_norm, std_Fy_error_norm, avg_Z_error_norm, std_Z_error_norm, a0s, a1s, b0s, b1s)
%Print results tables for copying into LaTeX
params = all_params(1);
fprintf('\n')
forcetable = [betas', avg_Fx_error_norm', std_Fx_error_norm', avg_Fy_error_norm', std_Fy_error_norm', avg_Z_error_norm', std_Z_error_norm'];
matrix2latex(forcetable, 1);
fprintf('\n')
abtable = [betas', a0s(1,:)', a1s(1,:)', b0s(1,:)', b1s(1,:)', ones(size(betas))'*params.rover.zeta];
matrix2latex(abtable, 2);
fprintf('\n')
end

function [fig, ax, tiles] = force_sinkage_plot(all_params, betas, Vrys, slip_data, Fx_data, Fy_data, Z_data, Fx_approx, Fy_approx, h, sinkage_found, fig, ax, tiles)
%Plot model computed forces and sinkages vs measured data
params = all_params(1);
if params.options.trench_on == 1 && params.options.ishigami_sidewall == 0
    linecolor = cmuColor('red-web'); %red plot for model with soil flow
    linestyle = ':';
elseif params.options.trench_on == 0 && params.options.ishigami_sidewall == 0
    linecolor = cmuColor('palladian-green'); %green plot for model without soil flow
    linestyle = '-';
elseif params.options.trench_on == 1 && params.options.ishigami_sidewall == 1
    linecolor = cmuColor('sky-blue'); %blue plot for Ishigami version
    linestyle = '--';
else
    linecolor = cmuColor('black');
    linestyle = '.-';
    fprintf('Check that options are set correctly')
end
dotcolor = cmuColor('dark-gray');

[~,N] = size(betas);
N = 5; %Override to exclude beta = 75, 90

if exist('fig') && ~isgraphics(fig)
    clear fig;
elseif exist('fig') && isgraphics(fig)
    fig = figure(1);
end
if ~exist('fig')
    fig = figure();
    fig.Position = [200 100 1100 600];
    tiles = tiledlayout(3, N, 'TileSpacing', 'tight');
    tiles.Padding = 'compact';
    ax = gobjects(3,N);
    if params.rover.h_g == 0
        title(tiles, 'Smooth Wheel Forces and Sinkages by Slip Angle', 'Interpreter','latex', 'FontSize', 20)
    else
        title(tiles, 'Grousered Wheel Forces and Sinkages by Slip Angle', 'Interpreter','latex', 'FontSize', 20)
    end
    xlabel(tiles, 'Slip Ratio', 'Interpreter','latex', 'FontSize', 16)
    for i=1:N
        ax(1,i) = nexttile;
        title([num2str(betas(i)), '$^\circ$'])
        axis([-1.1 1 -50 25])
        hold on
        ylabel('$F_x$ (N)', 'FontSize', 16)
        if i>1
            set(gca, 'YColor', 'none')
        end
        set(gca, 'XColor', 'none')
    end
    for i=1:N
        ax(2,i) = nexttile;
        axis([-1.1 1 -10 50])
        hold on
        ylabel('$F_y$ (N)', 'FontSize', 16)
        if i>1
            set(gca, 'YColor', 'none')
        end
        set(gca, 'XColor', 'none')
    end

    for i=1:N
        ax(3,i) = nexttile;
        axis([-1.1 1 0 90])
        hold on
        ylabel('Sinkage (mm)', 'FontSize', 16)
        if i>1
            set(gca, 'YColor', 'none')
        end
    end
end

%Plot asterisk for points where we failed to find the sinkage
failed_sinkages = 0;
for j = 1:length(Vrys)
    for k=1:N
        i = (j-1)*N + k;
        params = all_params(i);
        if sinkage_found(j,k) == 0 && params.options.ishigami_sidewall == 0
            marker = '*k';
            failed_sinkages = failed_sinkages + 1;
            plot(ax(1,k), slip_data(j,k), 1000*Fx_approx(j,k), marker, 'MarkerEdgeColor', dotcolor)

            plot(ax(2,k), slip_data(j,k), 1000*Fy_approx(j,k), marker, 'MarkerEdgeColor', dotcolor)

            plot(ax(3,k), slip_data(j,k), 1000*h(j,k), marker, 'MarkerEdgeColor', dotcolor)
        end
    end
end
%Plot the calculated forces as lines
for i=1:N
    % Plot x force model
    if params.options.ishigami_sidewall == 0
        plot(ax(1,i), slip_data(:,i), 1000*Fx_approx(:,i), 'LineStyle', linestyle, 'Color', linecolor, 'LineWidth', 2)
    end
end

for i=1:N
    %plot y force model
    plot(ax(2,i), slip_data(:,i), 1000*Fy_approx(:,i), 'LineStyle', linestyle, 'Color', linecolor, 'LineWidth', 2)
end

for i=1:N
    %plot z force model
    if params.options.ishigami_sidewall == 0
        plot(ax(3,i), slip_data(:,i), 1000*h(:,i), 'LineStyle', linestyle, 'Color', linecolor, 'LineWidth', 2)
    end
end

%Plot the raw data as dots
for i=1:N
    % Plot x force data
    plot(ax(1,i), slip_data(:,i), Fx_data(:,i), 'o', 'MarkerEdgeColor', dotcolor)
end

for i=1:N
    %plot y force data
    plot(ax(2,i), slip_data(:,i), Fy_data(:,i), 'o', 'MarkerEdgeColor', dotcolor)
end

for i=1:N
    %plot z force data
    plot(ax(3,i), slip_data(:,i), Z_data(:,i), 'o', 'MarkerEdgeColor', dotcolor)
end

bigax = axes(fig, 'visible', 'off');
xlabel(bigax, 'Slip Ratio')
hold off
end

function legend_plot()
%Generate legend for forces & sinkages plot
legf = figure();
legf.Position = [200 100 1100 50];
leg = zeros(3, 1);
hold on
plot(NaN, NaN, 'LineStyle', ':', 'Color', cmuColor('red-web'), 'LineWidth', 2);
plot(NaN, NaN, 'LineStyle', '-', 'Color', cmuColor('palladian-green'), 'LineWidth', 2);
plot(NaN, NaN, 'LineStyle', '--', 'Color', cmuColor('sky-blue'), 'LineWidth', 2);
leg = legend('With soil flow','Without soil flow','Ishigami 2007');
leg.Location = 'north';
leg.Orientation = 'horizontal';
axis off
hold off
legend 'boxoff'
end

function tuning_params_plot(all_params, betas, a0s, a1s)
%Plot a0 and a1 tuning parameters as a function of beta
if length(betas) > 1
    a0s_fit = fit(betas', a0s(1,:)', 'poly1');
    a1s_fit = fit(betas', a1s(1,:)', 'poly1');
end

f= figure();
f.Position = [200 200 600 350];
params = all_params(1);
tiles = tiledlayout(1,2);
if params.rover.h_g == 0
    title(tiles, 'Smooth Wheel Tuned Parameters as a Function of Slip Angle', 'Interpreter','latex')
else
    title(tiles, 'Grousered Wheel Tuned Parameters as a Function of Slip Angle', 'Interpreter','latex')
end

nexttile
coeffs = coeffvalues(a0s_fit);
hold on
title('$a_0$')
plot(betas, a0s(1,:), 'Color', cmuColor('dark-gray'))
plot(betas, coeffs(2) + coeffs(1)*betas, 'Color', cmuColor('red-web'))
hold off
xticks([0 45 90])
xlabel('$\beta$ (degrees)', 'FontSize', 14)
ylabel('$a_0$', 'FontSize', 14)
axis([0 90 -.1 1])

nexttile
coeffs = coeffvalues(a1s_fit);
hold on
title('$a_1$')
plot(betas, a1s(1,:), 'Color', cmuColor('dark-gray'))
plot(betas, coeffs(2) + coeffs(1)*betas, 'Color', cmuColor('red-web'))
hold off
xticks([0 45 90])
xlabel('$\beta$ (degrees)', 'FontSize', 14)
ylabel('$a_1$', 'FontSize', 14)
axis([0 90 -.1 1])

leg = legend('Tuned values', 'Linear fit') ;
leg.Layout.Tile ='south';
leg.Orientation = 'horizontal';
end