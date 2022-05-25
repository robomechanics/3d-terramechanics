%Plots all measured data points
load(['trench_off_grouser.mat']); %load data
grousered_results = tuned_results;

load('trench_off_smooth.mat');
smooth_results = tuned_results;

betas = [0 15 30 45 60 75 90];
Vrys = [0 3 5 8 10 12 20 40 100];

[g_Fx_data, g_Fy_data, g_Fz_data, g_Z_data, g_slip_data, g_V_data, g_Vry_data, g_angle_data, g_a0s, g_a1s, g_b0s, g_b1s, g_theta_m0s, g_theta_r0s] = average_data_to_model(grousered_results, betas, Vrys);
[s_Fx_data, s_Fy_data, s_Fz_data, s_Z_data, s_slip_data, s_V_data, s_Vry_data, s_angle_data, s_a0s, s_a1s, s_b0s, s_b1s, s_theta_m0s, s_theta_r0s] = average_data_to_model(smooth_results, betas, Vrys);

forcefig = figure();
sgtitle('Wheel Forces (from Z force balance)')
subplot(1,4,1)
title('Fx')
axis([-1.1 1 -50 25])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,4,2)
title('Fy')
axis([-1.1 1 -80 10])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,4,3)
title('Fz (Load)')
axis([-1.1 1 -25 40])
hold on
xlabel('Slip Ratio')
ylabel('Force (N)')

subplot(1,4,4)
title('Sinkage')
xlabel('Slip')
ylabel('Depth (mm)')
hold on

colors = [cmuColor('red-web'); cmuColor('gold'); cmuColor('teal'); cmuColor('sky-blue'); cmuColor('palladian-green'); cmuColor('blue-thread'); cmuColor('scots-rose')];

%Plot the raw data as dots
subplot(1,4,1) % Plot x force data
for i=1:length(betas)
    plot(g_slip_data(:,i), g_Fx_data(:,i), '^', 'Color', colors(i,:))
    plot(s_slip_data(:,i), s_Fx_data(:,i), 'o', 'Color', colors(i,:))
end

subplot(1,4,2) %plot y force data
for i=1:length(betas)
    plot(g_slip_data(:,i), g_Fy_data(:,i), '^', 'Color', colors(i,:))
    plot(s_slip_data(:,i), s_Fy_data(:,i), 'o', 'Color', colors(i,:))
end

subplot(1,4,3) %plot z force data
for i=1:length(betas)
    plot(g_slip_data(:,i), g_Fz_data(:,i), '^', 'Color', colors(i,:))
    plot(s_slip_data(:,i), s_Fz_data(:,i), 'o', 'Color', colors(i,:))
end

subplot(1,4,4); %plot sinkage data
for i=1:length(betas)
    plot(g_slip_data(:,i), g_Z_data(:,i), '^', 'Color', colors(i,:))
    plot(s_slip_data(:,i), s_Z_data(:,i), 'o', 'Color', colors(i,:))
end

leg(1) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('red-web'));
leg(2) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('gold'));
leg(3) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('teal'));
leg(4) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('sky-blue'));
leg(5) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('palladian-green'));
leg(6) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('blue-thread'));
leg(7) = plot(NaN, NaN, 'o', 'MarkerEdgeColor', cmuColor('scots-rose')); 
leg(8) = plot(NaN, NaN, 'ok');
leg(9) = plot(NaN, NaN, '^k');
legg = legend(leg, '0', '15', '30', '45', '60', '75', '90', 'Smooth', 'Grousered', 'Location', 'East');
title(legg, 'Slip Angle');
hold off

