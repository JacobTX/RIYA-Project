function F = disc_spring_force(x,ht_ratio,tau)

    E = 200*10^9/1e6;  % N/mm^2
    a = 34.5/2;        % Outer diameter, mm
    b = 22.4/2;        % Inner diameter, mm
    h = ht_ratio * tau; % height, mm

    alpha = a/b; % alpha (ratio of outer to inner diameter)
    M = ((alpha+1)/(alpha-1)-2/log(alpha))*tau; % M 
    N = (tau^3/6)*log(alpha); % N
    
    F = (E*pi/a^2) * (alpha/(alpha-1))^2 * x * ((h - x)*(h - x/2)*M + N); % Force exerted by disc-spring
end