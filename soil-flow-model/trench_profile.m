function [sp1,sp2,wp,d, p2, w, q, Atot, ip, ip_left, ip_right, int_profile] = trench_profile(params)

%get depths and widths of trench profile
[d, hl2, hr2, p2, hl1, hr1, w, q, m, d1, Atot] = trench_depth(params);

h = params.terr.h;
beta = params.state.beta;
r = params.rover.r;
b = params.rover.b;
phi = params.soil.phi;


%generate points for soil profile
if (p2 == 0) %pointy bottomed trench case
    if (q == 0) %no step
        soil_points1 = [-p2 - 4*r, 0;
            -w/2 - hl1/tan(phi), 0;
            -w/2 - hl1/tan(phi) + hl2/tan(phi), hl2;
            -w/2 - hl1/tan(phi) + (2*hl2 + d)/tan(phi), -d;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + hr2)/tan(phi), hr2;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + 2*hr2)/tan(phi), 0;
            p2 + 4*r, 0
            ];
        soil_points3 = soil_points1;
    else %step with pointy bottom 
        soil_points1 = [-p2 - 4*r, 0;
            -w/2 - hl1/tan(phi), 0;
            -w/2 - hl1/tan(phi) + hl2/tan(phi), hl2;
            -w/2 + (-hl1 + 2*hl2 + d - m)/tan(phi), -d + m;
            -w/2 + (-hl1 + 2*hl2 + d - m)/tan(phi) + q, -d + m;
            -w/2 + (-hl1 + 2*hl2 + d)/tan(phi) + q, -d;
            -w/2 + (-hl1 + 2*hl2 + 2*d + hr2)/tan(phi) + q, hr2;
            -w/2 + (-hl1 + 2*hl2 + 2*d + 2*hr2)/tan(phi) + q, 0;
            p2 + 4*r, 0];
        
        
        soil_points3 = [-p2 - 4*r, 0;
            -w/2 - hl1/tan(phi), 0;
            -w/2 - hl1/tan(phi) + hl2/tan(phi), hl2;
            -w/2 - hl1/tan(phi) + (2*hl2 + d1)/tan(phi), -d1;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d1 + hr2)/tan(phi), hr2;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d1 + 2*hr2)/tan(phi), 0;
            p2 + 4*r, 0
            ];
    end
else %flat bottomed trench case
    if (q == 0) %no step 
        soil_points1 = [-p2 - 3*r, 0;
            -w/2 - hl1/tan(phi), 0;
            -w/2 - hl1/tan(phi) + hl2/tan(phi), hl2;
            -w/2 - hl1/tan(phi) + (2*hl2 + d)/tan(phi), -d;
            -w/2 - hl1/tan(phi) + (2*hl2 + d)/tan(phi) + p2, -d;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + hr2)/tan(phi) + p2, hr2;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + 2*hr2)/tan(phi) + p2, 0;
            p2 + 4*r, 0];
        
        soil_points3 = soil_points1;
        
    else %flat bottom with step
        soil_points1 = [-p2 - 3*r, 0;
            -w/2 - hl1/tan(phi), 0;
            -w/2 - hl1/tan(phi) + hl2/tan(phi), hl2;
            -w/2 - hl1/tan(phi) + (2*hl2 + d - m)/tan(phi), -d + m;
            -w/2 - hl1/tan(phi) + (2*hl2 + d - m)/tan(phi) + q, -d + m;
            -w/2 - hl1/tan(phi) + (2*hl2 + d)/tan(phi) + q, -d;
            -w/2 - hl1/tan(phi) + (2*hl2 + d)/tan(phi) + p2, -d;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + hr2)/tan(phi) + p2, hr2;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + 2*hr2)/tan(phi) + p2, 0;
            p2 + 4*r, 0];
        
        soil_points3 = [-p2 - 3*r, 0;
            -w/2 - hl1/tan(phi), 0;
            -w/2 - hl1/tan(phi) + hl2/tan(phi), hl2;
            -w/2 - hl1/tan(phi) + (2*hl2 + d)/tan(phi), -d;
            -w/2 - hl1/tan(phi) + (2*hl2 + d)/tan(phi) + p2, -d;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + hr2)/tan(phi) + p2, hr2;
            -w/2 - hl1/tan(phi) + (2*hl2 + 2*d + 2*hr2)/tan(phi) + p2, 0;
            p2 + 4*r, 0];
    end
end

%generate points describing wheel profile

abeta = mod(beta, pi);
p1 = b*cos(beta);
t = linspace(0,pi,15)';
w_right = [p1/2 + r*sin(beta)*sin(t), r-h + r*cos(t)];
w_left  = [-1*w_right(:,1), w_right(:,2)];
w_back  = [flipud(w_left(:,1)) + b*cos(beta), flipud(w_left(:,2));
            w_right(:,1) - b*cos(beta), w_right(:,2)];
wheel_points = [-p1/2, -h;
    flipud(w_right);
    (w_left); (w_back)];
%w_center = (w_left + flipud(w_back))./2;
%w_center = w_center(1:end-1,:o);
w_left_center = [-p1/2, r-h];
w_right_center = [p1/2, r-h];

%Generate points describing intermediate soil profile
soil_points2 = [-p1 - 3*r, 0;
    -w/2 - hl1/tan(phi), 0;
    -w/2, hl1;
    -w/2, 0;
    0, 0;
    w/2, 0;
    w/2, hr1;
    w/2 + hr1/tan(phi), 0;
    w + 3*r, 0];

%flip profile is beta is negative
if beta<0
    soil_points1(:,1) = -soil_points1(:,1);
    soil_points2(:,1) = -soil_points2(:,1);
    soil_points3(:,1) = -soil_points3(:,1);
    soil_points1 = flipud(soil_points1);
    soil_points2 = flipud(soil_points2);
    soil_points3 = flipud(soil_points3);
end

%Find intersection of rear of wheel with soil
ip_left = []; %list of all intersection points
for i = 1:length(soil_points1) - 1
    line = soil_points1(i:i+1,:); %TODO make general
    new_ip = intersect_oval(line, w_left_center, r, beta);
    [n, ~] = size(new_ip);
    temp = [];
    if n > 0
        for j=1:n
            if sign(new_ip(j,1) - w_left_center(1)) == sign(-beta) %only take intersection points on the rear of the wheel
                temp = [temp; new_ip(j,:)]; 
            elseif (abs((new_ip(j,1) - w_left_center(1)))) < 10e-9 %Catch numerical error when the intersection point is right below wheel center, always true for beta=0 and some flat bottom trenches
                temp = [temp; new_ip(j,:)];
            end
        end
    end
    ip_left = [ip_left; temp];
end 
%if there are multiple intersection points, take the one highest up on the
%wheel
if isempty(ip_left)
    fprintf("error: exit angle not found")
else
    [~, max_ind] = max(ip_left(:,2));
    ip_left = ip_left(max_ind,:);
end



ip_right = []; %list of all intersection points
for i = 1:length(soil_points1) - 1
    line = soil_points1(i:i+1,:); %TODO make general
    new_ip = intersect_oval(line, w_right_center, r, beta);
        [n, ~] = size(new_ip);
    temp = [];
    if n > 0
        for j=1:n
            if (sign(new_ip(j,1) - w_right_center(1)) == sign(-beta)) 
                temp = [temp; new_ip(j,:)];
            elseif (abs((new_ip(j,1) - w_right_center(1)))) < 10e-9  %Catch numerical error when the intersection point is right below wheel center, always true for beta=0 and some flat bottom trenches
                temp = [temp; new_ip(j,:)];
            end
        end
    end
    ip_right = [ip_right; temp];
end  
if isempty(ip_right)
    fprintf("error: exit angle not found\n")
else
    [~, max_ind] = max(ip_right(:,2));
    ip_right = ip_right(max_ind,:);
end

%error handling for when the wheel is completely submerged
if isempty(ip_right)
    if isempty(ip_left)
        ip_left = wheel_points(17,:);
        ip_right = wheel_points(16,:);
    else
        ip_right = ip_left;
    end
elseif isempty(ip_left)
    ip_left = ip_right;
end


%Take left and right intersection points and find soil profile between them
int_profile = ip_left;
for i=1:length(soil_points1)
    if soil_points1(i,1) >= ip_left(1) && soil_points1(i,1) <= ip_right(1)
        int_profile = [int_profile; soil_points1(i,:)];
    end
end
int_profile = [int_profile; ip_right];

%Find weighted average intersection point
A = trapz(int_profile(:,1), int_profile(:,2));
ip = [(ip_left(1) + ip_right(1))/2, A/(ip_right(1) - ip_left(1))]; 

%return points for plotting
sp1 = soil_points1;
sp2 = soil_points2;
wp = wheel_points;
sp3 = soil_points3;



%sp2 = sp3; %Uncomment to display intermediate profile
end



