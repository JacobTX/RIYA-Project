%% Function used to calculate the spring force (F) for a given deflection (x)

% Based on non-linear analytical model
function F = disc_spring_force(x,ht_ratio,tau)
    k = 1.05;
    E = k * 200*10^9/1e6;  % Young's Modulus (N/mm^2)
    a = 34.5/2;        % Outer diameter (mm)
    b = 22.4/2;        % Inner diameter (mm)
    h = ht_ratio * tau; % Height (mm)

    alpha = a/b; % alpha (ratio of outer to inner diameter)
    M = ((alpha+1)/(alpha-1)-2/log(alpha))*tau; % M 
    N = (tau^3/6)*log(alpha); % N
    
    F = (E*pi/a^2) * (alpha/(alpha-1))^2 * x * ((h - x)*(h - x/2)*M + N); % Force exerted by disc-spring in Newton
end