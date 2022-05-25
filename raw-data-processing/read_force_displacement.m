function [k, n, Zs, Pzs, Fzs, expfit] = read_force_displacement(foldername, filename, min_time, max_time, disp_plot, DST, z_offset, Fz_offset)
if ~isempty(Fz_offset)
    Fz_zero = Fz_offset;
else
    Fz_zero = 0;
end

zerobag = rosbag(sprintf('%s%s_zero.bag', foldername, filename));

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
zavg(3) = mean(zeroFT(:,3)) - Fz_zero;
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




[X, Vx, Vry, angle, Z, sys_time, arduino_time] = read_speed_data(speedfilename, 0, 10000, 0);
if (DST == 0)
    sys_time = sys_time + ones(size(sys_time))*(5*3600 - bag.StartTime); %Correct for time zone error and shift time to start at the same zero as for the bag file non-DST
else
    sys_time = sys_time + ones(size(sys_time))*(4*3600 - bag.StartTime); %Correct for time zone error and shift time to start at the same zero as for the bag file DST
end

if ~isempty(z_offset)
    Z_zero = z_offset;
else
    Z_zero = 0; %sinkage of wheel at starting position
end

Z = Z - Z_zero;

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


%Plot data including lines to show the extent of the cropping
if disp_plot  
    figure
    hold on
    plot(time, movmean(RMF(:,1:3), 10));
    plot(sys_time, Z);
    xline(sys_time(sys_lower_ind));
    xline(sys_time(sys_upper_ind));
    legend('Fx', 'Fy', 'Fz', 'Z', 'Location', 'Southwest')
    ylim([-150 150])
    hold off
end


%Now crop the data

RMF_crop = RMF(ft_lower_ind:ft_upper_ind,:);
time_crop = time(ft_lower_ind:ft_upper_ind);
X_crop = X(sys_lower_ind:sys_upper_ind);
Vx_crop = Vx(sys_lower_ind:sys_upper_ind);
Vry_crop = Vry(sys_lower_ind:sys_upper_ind);
angle_crop = angle(sys_lower_ind:sys_upper_ind);
Z_crop = Z(sys_lower_ind:sys_upper_ind);
sys_time_crop = sys_time(sys_lower_ind:sys_upper_ind);
arduino_time_crop = arduino_time(sys_lower_ind:sys_upper_ind);

%Interpolate the force and displacement data to be the same size
Fz_crop = -movmean(RMF_crop(:,3), 10)/1000; %flip sign and switch from N to kN TODO check units
Z_crop_interp = -interp1q(sys_time_crop, Z_crop, time_crop)/1000; %interpolate, flip sign and switch from mm to m TODO check units
Z_crop_interp = Z_crop_interp(10:end-2); %crop off NaN values at the ends
Fz_crop = Fz_crop(10:end-2); %crop length of Fzs to match 

deflection_compensation = 0;
if deflection_compensation
    linfit = testbed_deflection_analysis(0); %get compensation for deflection of the testbed
    Z_crop_interp = Z_crop_interp - linfit(Fz_crop); %Compensate for testbed deflection
end
    
%Convert forces to pressures 
r = .0580; %hardocde value for r
b = .050; %hardcode value for b


theta = acos((r*ones(length(Z_crop_interp),1)-Z_crop_interp)./(r*ones(length(Z_crop_interp),1)));
A = b*2*r*sin(theta);
Pz_crop = Fz_crop./A;

%Construct a pressure-displacement curve from the cropped data
expfit = fit(double(Z_crop_interp), double(Pz_crop), 'power1');

if disp_plot
    figure()
    hold on
    plot(Z_crop_interp, Pz_crop)
    plot(expfit, Z_crop_interp, Pz_crop)
    xlabel('Sinkage(m)')
    ylabel('Pressure(kPa)')
    hold off
end

coeffs = coeffvalues(expfit);
k = coeffs(1);
n = coeffs(2);

Zs = Z_crop_interp;
Pzs = Pz_crop;
Fzs = Fz_crop;


