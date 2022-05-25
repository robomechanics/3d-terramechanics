function [X, Vx, Vry, angle, Z, sys_time, arduino_time, Ry] = read_speed_data(filename, min_x, max_x, disp_plot)
a = readmatrix(filename);
[m,~] = size(a);
lower_ind = 1;
upper_ind = m;
k=1;

while a(k,1) <= max_x && k < m
    if a(k,1) <= min_x
        lower_ind = k+1;
    end
    k = k+1;
end
upper_ind = k;

a_cropped = a(lower_ind:upper_ind,:);

X = a_cropped(:,1);
Vx = a_cropped(:,2);
Vry = a_cropped(:,3);
angle = a_cropped(:,4);
Z = a_cropped(:,5);
arduino_time = (a_cropped(:,6)-ones(upper_ind - lower_ind + 1,1)*a_cropped(1,6))/1000; 
sys_time = a_cropped(:,7);
[~,n] = size(a);
if n == 8
    Ry = a_cropped(:,8); %Needed for backwards compatibility with old data format
else
    Ry = []; %Needed for backwards compatibility with old data format
end

if disp_plot
    figure()
    hold on
    plot(sys_time - sys_time(1), Vx, 'DisplayName', 'Vx');
    plot(sys_time - sys_time(1), Vry, 'DisplayName', 'Vry');
    plot(sys_time - sys_time(1), angle, 'DisplayName', 'Angle');
    plot(sys_time - sys_time(1), Z, 'DisplayName', 'Z');
    if n==8
        plot(sys_time - sys_time(1), Ry, 'DisplayName', 'Ry');
    end
    hold off
    legend
end
    
