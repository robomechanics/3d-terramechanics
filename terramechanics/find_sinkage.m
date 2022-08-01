%Finds the sinkage for a given state by performin a binary search over 
%sinkage to find a sinkage that balances the applied vertical load Fz
function [hz, params] = find_sinkage(params)

W = params.rover.W;
r = params.rover.r;
res = .0001; %desired resolution for h in m
h_min = res;
h_max = 1.5*r; %Max is set to 1.5r instead of 2r to avoid nonmonotonically increasing section of Fz above 1.6r. 
Fz = 0;
err = .001*W;
params.terr.sinkage_found = 1; %assume we will be able to find the sinkage until proven otherwise

%Use binary search to find sinkage to desired resolution
while abs(W-Fz) > abs(err) %.1% error on force balance
    %check to see if force balance has been met
    h = (h_max + h_min)/2; %update the value to test
    params.terr.h = h;
    params = update_thetas(params);
    [~, ~, Fz, params] = forces(params);
    if imag(Fz) ~= 0 %Generally corresponds to theta_m0 > theta_f, so we need a larger sinkage
        h_min = h;
    elseif Fz < W
        h_min = h;
    else
        h_max = h;
    end
    if h_max - h_min < .00000001
        v_x = params.state.v_x;
        v_y = params.state.v_y;
        V = sqrt(v_x^2 + v_y^2);
        if params.state.w*params.rover.r_s >= V %slip
            s = (params.state.w*params.rover.r_s - V)/(params.state.w*params.rover.r_s); 
        else %skid
            s = (params.state.w*params.rover.r_s - V)/(V); 
        end
        fprintf('Unable to find sinkage for beta = %.f s = %.1f, current Fz = %.4f, desired Fz = %.4f\n', rad2deg(params.state.beta), s, Fz, W)
        params.terr.sinkage_found = 0; %Note that we did not find the sinkage correctly here
        break
    end
end
params.terr.h = h;
hz = h;
    
    