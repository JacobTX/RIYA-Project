%% Function that represents the non-linear dynamics of the 2-spring stack

function dx = dyn_dspring_stack(t,x,m1,g,ht1_ratio,tau1,ht2_ratio,tau2,A,f_excit)

    % Base displacement
    amp_ptp = A; % Peak-to-Peak amplitude (mm)
    x_base_max = amp_ptp/2; % Amplitude of displacement mm
    f_base = f_excit; % Frequency of base excitation (Hz)
    w_base = 2 * pi * f_base; % Angular frequency of base excitation (rad/s)
    x_base = x_base_max * sin(w_base * t); % Base displacement (mm)
    xd_base = x_base_max * w_base * cos(w_base * t); % Base velocity (mm/s)

    % Driving Force
    F0 = 0; % in N
    w_drive = 10;
    F_drive = F0 * sin(w_drive * t);

    % Spring Forces
    fun = @(y)equality(y, x(1,1), x_base, ht1_ratio, tau1, ht2_ratio, tau2); % Spring 1 Force = Spring 2 Force
    y = fsolve(fun, 0); % Solve for displacement of spring 2
    % Tolerance in fsovle may be decreased if necessary
    
    F1_spring = disc_spring_force(x(1,1) - y, ht1_ratio, tau1); % Value of Spring 1 Force
    
    % Spring Forces
    %F1_spring = disc_spring_force_exp(x(1,1)); % Value of Spring 1 Force based on experimental model

    % ****** CHANGE DAMPING PARAMETERS HERE ******
   
    % Linear Viscous Damping Force
    c = 0.2;
    F_damp = c * (x(2,1) - xd_base);

    % Coulomb Friction
    mu = 0.03;
    F_stat = mu * (m1 * g) * tanh((x(2,1)-xd_base));

    % State evolution
    dx(1,1) = x(2,1); % Velocity in mm/s
    dx(2,1) = 1000 * (g + (F_drive - F1_spring - F_damp - F_stat)/m1); % Acceleration in mm/s^2
    dx(3,1) = xd_base; % Base excitation velocity in mm/s
    dx(4,1) = -0.001 * m1 * g * x(2,1)  + 0.001 * F1_spring * (x(2,1) - xd_base); % Potential energy of the system (Gravity + Spring)
end