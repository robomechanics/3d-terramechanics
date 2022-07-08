%Takes individual raw data points and groups by slip angle and slip ratio
%to calculate mean and standard deviation values for each (slip angle, slip
%ratio) pair. Also extracts the tuned constants and returns them in the
%same format.
function [Fxs, Fys, Fzs, Zs, slips, V_means, Vry_means, angle_means, a0s, a1s, b0s, b1s, theta_m0s, theta_r0s] = average_data_to_model(all_results, betas, Vrys)

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
V_means = zeros(length(betas), length(Vrys));
V_stds = zeros(length(betas), length(Vrys));
Vry_means = zeros(length(betas), length(Vrys));
Vry_stds = zeros(length(betas), length(Vrys));
angle_means = zeros(length(betas), length(Vrys));
angle_stds = zeros(length(betas), length(Vrys));
a0s = zeros(length(betas), length(Vrys));
a1s = zeros(length(betas), length(Vrys));
b0s = zeros(length(betas), length(Vrys));
b1s = zeros(length(betas), length(Vrys));
theta_m0s = zeros(length(betas), length(Vrys));
theta_r0s = zeros(length(betas), length(Vrys));

temp = [];
Fx = [];
Fy = [];
Fz = [];
Z = [];
slip = [];
V = [];
Vry = [];
angle = [];

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
            Fy = [Fy, datapt.avg_Fx];
            Fz = [Fz, -datapt.avg_Fz];
            Z = [Z, -datapt.avg_Z];
            slip = [slip, datapt.slip];
            V = [V, datapt.avg_Vx];
            Vry = [Vry, datapt.avg_Vry];
            angle = [angle, datapt.avg_angle];
            
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
        V_means(j,k) = mean(V);
        V_stds(j,k) = std(V);
        Vry_means(j,k) = mean(Vry);
        Vry_stds(j,k) = std(Vry);
        angle_means(j,k) = mean(angle);
        angle_stds(j,k) = std(angle);
        
        if isfield(datapt, 'a0') %Assume that if we have one of these variables tuned they are all present
            a0s(j,k) = datapt.a0;
            a1s(j,k) = datapt.a1;
            b0s(j,k) = datapt.b0;
            b1s(j,k) = datapt.b1;
            theta_m0s(j,k) = datapt.theta_m0;
            theta_r0s(j,k) = datapt.theta_r0;
        end
        
        temp = [];
        Fx = [];
        Fy = [];
        Fz = [];
        Z = [];
        slip = [];
        V = [];
        Vry = [];
        angle = [];
    end
end
    


Fxs = Fx_means';
Fys = Fy_means';
Fzs = Fz_means';
Zs = Z_means';
slips = slip_means';
V_means = V_means';
Vry_means = Vry_means';
angle_means = angle_means';
a0s = a0s';
a1s = a1s';
b0s = b0s';
b1s = b1s';
theta_m0s = theta_m0s';
theta_r0s = theta_r0s';
