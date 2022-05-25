%takes a line and checks to see if it intersects an oval defined by its
%center [x0, y0], radius r, and skew angle beta (such that the oval has a
%height of 2r and width of 2r*sin(beta)
function ip = intersect_oval(line, center, r, beta)
ip = [];

x0 = center(1);
y0 = center(2);
x1 = line(1,1);
y1 = line(1,2);
x2 = line(2,1);
y2 = line(2,2);

%Find intersections (note: not sure, but this might only return answers
%below r in y
a = (x2-x1)^2 + (y2-y1)^2*sin(beta)^2;
b = 2*(x2-x1)*(x1-x0) + 2*(y2-y1)*(y1-y0)*sin(beta)^2;
c = (x1-x0)^2 + (y1-y0)^2*sin(beta)^2 - r^2*sin(beta)^2;
tplus =  (-b + sqrt(b^2 - 4*a*c))/(2*a);
tminus = (-b - sqrt(b^2 - 4*a*c))/(2*a);
if tplus <= 1 && tplus >= 0
    xint = x1 + tplus*(x2-x1);
    yint = y1 + tplus*(y2-y1);
    ip = [ip; [xint yint]];
else
    ip = [ip;[]];
end
if tminus <= 1 && tminus >= 0
    xint = x1 + tminus*(x2-x1);
    yint = y1 + tminus*(y2-y1);
    ip = [ip; [xint yint]];
else
    ip = [ip;[]];
end
if ~isreal(ip)
    if abs(imag(ip))< 10e-6
        ip = real(ip);
    else
        ip = [];
    end
end

end