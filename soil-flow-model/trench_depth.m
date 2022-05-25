%see quals paper/presentation for notation

function [d, HL2, HR2, p2, HL1, HR1, w, q, m, d1, Atot] = trench_depth(params)
%%Define key geometry & driving parameters
r = params.rover.r; %radius of wheel
beta = params.state.beta; %slip angle
abeta = mod(abs(beta), pi); % maps angles to 0 to 90, assumes angle is always CW from North (trench_profile will take care of which side is leading edge)
if beta > pi/2
    abeta = pi - abeta;
end
h = params.terr.h; %sinkage of wheel
b = params.rover.b; %width of wheel
phi = params.soil.phi; %angle of repose for soil
p1 = b*cos(abeta);
w = p1 + 2*r*sin(abeta)*sin(acos((r-h)/r));
omega = params.state.w;
mu = params.rover.mu;
hg = params.rover.h_g;
vx = params.state.v_x;
zeta = params.rover.zeta;


%%set initial size of key soil profile features to 0 (defines trench as
%being totally flat)
hl1 = 0; %Initial height of sand pile to left of trench
hr1 = 0; %Initial height of sand pile to right of trench
hl2 = 0; %Final height of sand pile to left of trench
hr2 = 0; %Final height of sand pile to right of trench
d = 0; %depth of deepest part of trench
p2 = 0; %width of flat bottom of trench (p2 = 0 means trench is pointy)
q = 0; %width of step in trench
m = 0; %height of step in trench
Arot = 0; %amount of soil transported by rotation
Al1 = 0; %area of soil in left pile before flow, and before rotated soil is accounted for
Ar1 = 0; %area of soil in right pile before flow, and before rotated soil is accounted for
Al2 = 0; %area of soil in left pile before flow, and after rotated soil is accounted for
Ar2 = 0; %area of soil in right pile before flow, and after rotated soil is accounted for

%%Find initial heights of sand piles to right and left
%A_el = ellipse_int(r, abeta, h); %numerically determined area of elliptical sections
A_el = (w-p1)/3*h; %quadratic approximation to elliptical sections

tol = deg2rad(10); %Angle for which small angle approx holds well
%For large or small piles, sand is distributed evenly
if abeta < tol || (abeta > (pi/2-tol) && abeta < (pi/2+tol)) || abeta > (pi-tol) || (2*A_el <= b*h*cos(abeta))
    Al1 = .5*h*p1+ A_el;
    Ar1 = Al1;
    Atot = Ar1 + Al1;
    %otherwise, sand splits along leading corner of wheel
    %Note that there is a discontinuity in soil flow distribution
else
    Atot = (h*p1 + 2*A_el); %Projected area of wheel onto vertical plane
    if abeta > 0
        Al1 = 2*A_el;
        Ar1 = Atot - Al1;
    else
        Ar1 = A_el + h*p1*sin(abeta);
        Al1 = Atot - Ar1;
    end
end


%Find portion of sand transported by the grousers
if (hg == 0)
    hg = params.rover.r_s-r; %TODO: better approximation for shearing radius
    mu = 1;
end
Arot = 2*pi*r*omega*mu*b*hg/vx*zeta; %amount of soil transported by rotating grousers, assuming less than max amount is moved


%Low angle:
%reallocate sand from left and right piles to the grouser
%portion
if (abeta < tol)
    Arotmax = Atot;
    if Arot > Arotmax
        Arot = Arotmax;
    end
    Al2 = Al1 - Arot/2; %take half the roated sand from left
    Ar2 = Ar1 - Arot/2; %and half from the right
end

%Intermediate Angle:
%reallocate sand from left and right piles to the grouser
%portion

if (tol <= abeta) && (2*sqrt(2*r*h - h^2)*sin(abeta) < b*cos(abeta)) %this may never occur for some wheel geometries, the second condition is the point when the front and rear corners of the wheel alighn
    Arotmax = Ar1;
    if Arot > Arotmax
        Arot = Arotmax;
    end
    Al2 = Al1; %nothing happens to left pile
    Ar2 = Ar1 - Arot; % rotated soil gets removed from right pile
end

%High Angle:
%reallocate sand from left and right piles to the grouser
%portion
if (2*sqrt(2*r*h - h^2)*sin(abeta) > b*cos(abeta))
    Arotmax = Ar1;
    if Arot > Arotmax
        Arot = Arotmax;
    end
    Ar2 = Ar1 - Arot; % rotated soil gets removed from right pile
    
    %Put the sand transported by rotation at high angle in the correct pile
    Al2 = Al1 + Arot;
    Arot = 0;
end

%all cases:
%find height of the soil piles on either side of the wheel
hl1 = sqrt(2*tan(phi)*Al2);
hr1 = sqrt(2*tan(phi)*Ar2);




%%Find depth of trench

%Assumes the max trench depth occurs along the same line as the soil
%division
%Is this a crap assumption? Probably! Only actual experiments will tell.
%solution to quadratic, see notes for the math (need to LaTeX the math for
%later)

if hr1 == 0 % we have a different solution if there is no soil in the pile on the right (or the left, but we have set this up so soil is always removed from the right)
    d1 = 1/3*(2*sqrt((hl1 + (w+p1)/2*tan(phi))^2 + 3*Arot) - (hl1 + (w+p1)/2*tan(phi)));
elseif hl1 == 0
    d1 = 1/3*(2*sqrt((hl1 + (w-p1)/2*tan(phi))^2 + 3*Arot) - (hl1 + (w-p1)/2*tan(phi)));    
else
    A1 = -2;
    B1 = -2*(hl1+hr1+w*tan(phi));
    C1 = (hl1 + (w+p1)/2*tan(phi))^2 + (hr1 + (w-p1)/2*tan(phi))^2 + 4*Arot*tan(phi);%TODO double check math in the mornin, switched from -Arot to +Arot and it seemed to fix the issues with initially-pointy stepped trenches
    
    d1 = (-B1 - sqrt(B1^2-4*A1*C1))/(2*A1);
    d2 = (-B1 + sqrt(B1^2-4*A1*C1))/(2*A1); %this is unused
    
    if d1 < .00001 %for numerical errors
        d1 = 0;
    end  
end




%Checks to see if bottom of trench is flat
if d1 <= h + .00 %If not, returns pointy ttrench profile
    d = d1;
    if (abeta < tol) %driving straight, so trench max depth is centered
        temp = (hl1 + w/2*tan(phi));
        hl2 = temp - sqrt(temp^2 - 1/2*((temp^2 - Arot*tan(phi)))); % My hand calcs have +Arot, but -Arot works while + does not. TODO double check my math in "straight trench with arot calcs"
        hr2 = hl2;
        d = hl1 - 2*hl2 + w/2*tan(phi);
        
    else
        if hr1 == 0
            hl2 = sqrt(d^2 + Arot*tan(phi));
            hr2 = 0;
        elseif hl1 == 0
            hr2 = sqrt(d^2 + Arot*tan(phi));
            hl2 = 0;
        else
            hl2 = hl1/2 - d/2 + (w+p1)*tan(phi)/4; %This should be good, Arot is accounted for when calculating d.
            hr2 = hr1/2 - d/2 + (w-p1)*tan(phi)/4;
            hr2^2/tan(phi) + hl2^2/tan(phi) - d^2/tan(phi) + Arot;
        end
        
    end
    
    
    
%If trench depth is deeper than the wheel's sinkage, we have a flat-bottomed
%trench, recalculate soil pile heights, doing conservation of area on the
%left and right halves independently.
else
    d = h;
    if hl1 > 0
        A2 = 1;
        B2 = 2*h;
        C2 = d^2/2 - d*hl1 - hl1^2/2 - (w-p1)/2*d/3*tan(phi);
        hl2 = (-B2 + sqrt(B2^2-4*A2*C2))/(2*A2);
    else
        hl2 = 0;
    end
    
    if hr1 > 0
        A3 = 1;
        B3 = 2*h;
        C3 = d^2/2 - d*hr1 - hr1^2/2 - (w-p1)/2*d/3*tan(phi); %should be 0 when hr1= 0 TODO
        hr2 = (-B3 + sqrt(B3^2-4*A3*C3))/(2*A3); %should be 0 when hr1 = 0 TODO
    else
        hr2 = 0;
    end
    
    p2 = (hl2^2 + hr2^2 - d^2 + Arot*tan(phi))/(d*tan(phi)); %theoretically should work? Includes Arot now.
    %p2 = w + (hl1+hr1 - 2*hl2-2*hr2-2*d)/tan(phi); %this is wrong when hr1 = 0 because the trench gets a little wider at that edge during flow
    (hl2^2 + hr2^2 - d^2)/tan(phi) + Arot - d*p2; 
end


%%Now deal with soil transport for low and intermediate angles

%For low beta: find new depth of flat bottom
if (abeta < tol)
    if (p2 == 0)
        d = d - sqrt(Arot*tan(phi));
        p2 = 2*sqrt(Arot/tan(phi));
    else
        % new depth for initially flat bottom (note that these equations
        % simplify to those above when p2 == 0)
        d = d + p2*tan(phi)/2 - sqrt((p2*tan(phi)/2)^2 + Arot*tan(phi));
        p2 = sqrt(p2^2 + 4*Arot/tan(phi));
    end
    
    %For intermediate beta angles with nonzero Arot: find the dimensions of a
    %step, or the new depth of a flat bottom
elseif ((abeta < deg2rad(90)-tol) && (omega > 0))
    q = p1/2 + p2/2; %minimum width of step
    if (p2 == 0) %initially pointy trench
        d1 = d;
        m = Arot/q + q*tan(phi)/4; %height of step for initially pointy trench (I checked this it looks good)
        mmax = min((2*r*sin(abeta)+b*cos(abeta) - p2)*tan(phi)/2,hl2 + d); %max height of step
        if (m > mmax) %m is the height ALL the way to the deepest part of the trench
            m = mmax;
            q = 2*m/tan(phi) - sqrt((2*m/tan(phi))^2 - 4*Arot/tan(phi)); %checked this it good
            if (q >= (2*m/tan(phi)))
                q = 0; %no step
                m = 0;
                d = d - sqrt(Arot*tan(phi)); %just a flat bottom as in the low angle case
                p2 = 2*sqrt(Arot/tan(phi));
            else
                d = d - q/2*tan(phi); %checked this it good
                m = m - q/2*tan(phi); %checked this it good
            end
        else
            if (q >= (2*m/tan(phi)))
                q = 0; %no step
                m = 0;
                d = d - sqrt(Arot*tan(phi)); %just a flat bottom as in the low angle case
                p2 = 2*sqrt(Arot/tan(phi));
            else
                d = d - q/2*tan(phi);
                m = m - q/2*tan(phi);
            end
        end
        
    else %initially flat trench
        m = Arot/q; %height of step for initially flat trench
        mmax = min((2*r*sin(abeta)+b*cos(abeta) - p2)*tan(phi)/2, hl2 + d); %max height of step
        if (m > mmax)
            m = mmax; %if height would exceed max height, set height to max
            q = Arot/mmax;% and then recalculate step width
        end
        if (q > p2)%if the step is so wide it runs into the other side, recalculate width again
            if m*(m/tan(phi)+p2) > Arot %If this is true the step runs into the other side
                q = ((2*p2 + 4*m/tan(phi)) - sqrt((2*p2 + 4*m/tan(phi))^2 - 4*(p2^2 + 4*Arot/tan(phi))))/2; %this is correct, I checked it
                d = d - (q-p2)/2*tan(phi);
                m = m - (q-p2)/2*tan(phi); %still correct here, volume is conserved.
                p2 = 0;
            else %step runs into other side, is flat
                d = d + p2*tan(phi)/2 - sqrt((p2*tan(phi)/2)^2 + Arot*tan(phi));
                p2 = sqrt(p2^2 + 4*Arot/tan(phi));
                q = 0; %because technically we don't have a step anymore
                m = 0;
            end
        end
    end
end

%deal with the case where no soil is moved
if (hr1 == 0) && (hl1 == 0)
    hl2 = 0;
    hr2 = 0;
    d = 0;
    p2 = w;
    d1 = 0;
end

HL2 = hl2;
HR2 = hr2;
HL1 = hl1;
HR1 = hr1;




