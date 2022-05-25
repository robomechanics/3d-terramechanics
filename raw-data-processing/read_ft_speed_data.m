function results = read_ft_speed_data(foldername, filename, min_time, max_time, disp_plot, DST, def_fit)
Z_zero = 31.5; %sinkage of wheel at starting position
zerobag = rosbag(sprintf('%s%s_zero.bag', foldername, filename));
lin_fit = 1; %Should we do a linear fit to the data or just take the average

%Coefficients of polynomial fit to Fx, Fy, Fz
polyx = [];
polyy = [];
polyz = [];

bag = rosbag(sprintf('%s%s.bag', foldername, filename));
speedfilename = sprintf('%s%s.txt', foldername, filename);

%Actual data
RMtopic = select(bag,'Topic','/rokubimini_cosmo/forcetorque_readings');
RMStruct = readMessages(RMtopic,'DataFormat','struct');

%Zero file
Ztopic = select(zerobag,'Topic','/rokubimini_cosmo/forcetorque_readings');
ZStruct = readMessages(Ztopic, 'DataFormat','struct');

% Create a an Array of all 6 measurements history
zeroFT(:,1) = cellfun(@(m) m.Wrench.Wrench.Force.X,ZStruct);
zeroFT(:,2) = cellfun(@(m) m.Wrench.Wrench.Force.Y,ZStruct);
zeroFT(:,3) = cellfun(@(m) m.Wrench.Wrench.Force.Z,ZStruct);
zeroFT(:,4) = cellfun(@(m) m.Wrench.Wrench.Torque.X,ZStruct);
zeroFT(:,5) = cellfun(@(m) m.Wrench.Wrench.Torque.Y,ZStruct);
zeroFT(:,6) = cellfun(@(m) m.Wrench.Wrench.Torque.Z,ZStruct);
zts = cellfun(@(m) m.Wrench.Header.Stamp,ZStruct);

zsec = [zts.Sec] - zerobag.StartTime;
ztime = (cast(zsec, 'single') + cast([zts.Nsec], 'single')/10^9);

zavg(1) = mean(zeroFT(:,1));
zavg(2) = mean(zeroFT(:,2));
zavg(3) = mean(zeroFT(:,3));
zavg(4) = mean(zeroFT(:,4));
zavg(5) = mean(zeroFT(:,5));
zavg(6) = mean(zeroFT(:,6));

% Create a an Array of all 6 measurements history
RMF(:,1) = cellfun(@(m) m.Wrench.Wrench.Force.X,RMStruct) - zavg(1);
RMF(:,2) = cellfun(@(m) m.Wrench.Wrench.Force.Y,RMStruct) - zavg(2);
RMF(:,3) = cellfun(@(m) m.Wrench.Wrench.Force.Z,RMStruct) - zavg(3);
RMF(:,4) = cellfun(@(m) m.Wrench.Wrench.Torque.X,RMStruct) - zavg(4);
RMF(:,5) = cellfun(@(m) m.Wrench.Wrench.Torque.Y,RMStruct) - zavg(5);
RMF(:,6) = cellfun(@(m) m.Wrench.Wrench.Torque.Z,RMStruct) - zavg(6);
ts = cellfun(@(m) m.Wrench.Header.Stamp,RMStruct);

sec = [ts.Sec] - bag.StartTime;
time = (cast(sec, 'single') + cast([ts.Nsec], 'single')/10^9)';




[X, Vx, Vry, angle, Z, sys_time, arduino_time, Ry] = read_speed_data(speedfilename, 0, 10000, 0);
if (DST == 0)
    sys_time = sys_time + ones(size(sys_time))*(5*3600 - bag.StartTime); %Correct for time zone error and shift time to start at the same zero as for the bag file non-DST
else
    sys_time = sys_time + ones(size(sys_time))*(4*3600 - bag.StartTime); %Correct for time zone error and shift time to start at the same zero as for the bag file DST
end
Z = Z - Z(1) + Z_zero;

%find indices for cropping data
%First for the force data
ft_lower_ind = 1;
ft_upper_ind = length(time);
for j=1:length(time)
    if time(j,1) <= min_time
        ft_lower_ind = j+1;
    end
    if time(j,1) <= max_time
        ft_upper_ind = j;
    end
end

%And then for the speed data
sys_lower_ind = 1;
sys_upper_ind = length(sys_time);
for k=1:length(sys_time)
    if sys_time(k) <= min_time
        sys_lower_ind = k+1;
    end
    if sys_time(k) <= max_time
        sys_upper_ind = k;
    end
end

%Crop the data
RMF_crop = RMF(ft_lower_ind:ft_upper_ind,:);
X_crop = X(sys_lower_ind:sys_upper_ind);
Vx_crop = Vx(sys_lower_ind:sys_upper_ind);
Vry_crop = Vry(sys_lower_ind:sys_upper_ind);
angle_crop = angle(sys_lower_ind:sys_upper_ind);
Z_crop = Z(sys_lower_ind:sys_upper_ind);
sys_time_crop = sys_time(sys_lower_ind:sys_upper_ind);
time_crop = time(ft_lower_ind:ft_upper_ind);
arduino_time_crop = arduino_time(sys_lower_ind:sys_upper_ind);

%Compensate for deflection of testbed in the Z direction
%Interpolate the force and displacement data to be the same size
Fz_interp = interp1q(time_crop, RMF_crop(:,3), sys_time_crop); %interpolate
Fz_interp = Fz_interp(2:end-1); %crop off NaN values
Z_crop = Z_crop(2:end-1); %crop size to match
Z_crop = Z_crop - def_fit(Fz_interp); %Compensate for testbed deflection

%Polynomial fit to Fx, Fy, Fz
if lin_fit
    polyx = polyfit(time(ft_lower_ind:ft_upper_ind), RMF(ft_lower_ind:ft_upper_ind,1), 1);
    polyy = polyfit(time(ft_lower_ind:ft_upper_ind), RMF(ft_lower_ind:ft_upper_ind,2), 1);
    polyz = polyfit(time(ft_lower_ind:ft_upper_ind), RMF(ft_lower_ind:ft_upper_ind,3), 1);
end

%Plot data including lines to show the extent of the cropping
if disp_plot  
    figure
    hold on
    plot(time, movmean(RMF(:,1:3), 10));
    plot(sys_time, Vx);
    plot(sys_time, Vry);
    plot(sys_time, Z);
    plot(sys_time_crop(2:end-1), Z_crop, '--');
    plot(sys_time, angle);
    if lin_fit 
        plot(time(ft_lower_ind:ft_upper_ind), polyval(polyx,time(ft_lower_ind:ft_upper_ind)), '-k')
        plot(time(ft_lower_ind:ft_upper_ind), polyval(polyy,time(ft_lower_ind:ft_upper_ind)), '-k')
        plot(time(ft_lower_ind:ft_upper_ind), polyval(polyz,time(ft_lower_ind:ft_upper_ind)), '-k') 
    end
    xline(sys_time(sys_lower_ind));
    xline(sys_time(sys_upper_ind));
    legend('Fx', 'Fy', 'Fz', 'Vx', 'Vry', 'Z', 'Adjusted Z', 'Angle', 'Location', 'Southwest')
    ylim([-150 150])
    hold off
end

%Now calculate the average forces, torques, slip, and slip angle for the
%cropped data

%Forces and torques:
results.avg_Fx = mean(RMF_crop(:,1));
results.std_Fx = std(RMF_crop(:,1));
results.avg_Fy = mean(RMF_crop(:,2));
results.std_Fy = std(RMF_crop(:,2));
results.avg_Fz = mean(RMF_crop(:,3));
results.std_Fz = std(RMF_crop(:,3));

%Speeds and positions:
results.avg_Vx = mean(Vx_crop);
results.std_Vx = std(Vx_crop);
results.avg_Vry = mean(Vry_crop);
results.std_Vry = std(Vry_crop);
results.avg_angle = mean(angle_crop);
results.std_angle = std(angle_crop);
results.avg_Z = mean(Z_crop);
results.std_Z = std(Z_crop);
results.polyx = polyx;
results.polyy = polyy;
results.polyz = polyz;


if(disp_plot)
    fprintf(sprintf('Mean: Fx: %.2f   Fy: %.2f    Fz: %.2f\n', results.avg_Fx, results.avg_Fy, results.avg_Fz))
    fprintf(sprintf('Std: Fx: %.2f   Fy: %.2f    Fz: %.2f\n', results.std_Fx, results.std_Fy, results.std_Fz))
    fprintf(sprintf('Mean: Vx: %.2f   Vry: %.2f    angle: %.2f  Z: %.2f\n', results.avg_Vx, results.avg_Vry, results.avg_angle, results.avg_Z))
    fprintf(sprintf('Std: Vx: %.2f   Vry: %.2f    angle: %.2f  Z: %.2f\n', results.std_Vx, results.std_Vry, results.std_angle, results.std_Z))
end
