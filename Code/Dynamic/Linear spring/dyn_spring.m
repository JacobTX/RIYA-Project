function dx = dyn_spring(t,x,m,g,b)

    % Base displacement
    x_base_max = 1;
    w_base = 15;
    x_base = x_base_max * sin(w_base * t);
    xd_base = x_base_max * w_base * cos(w_base * t);

    % Driving Force
    F0 = 0;
    w_drive = 10;
    F_drive = F0 * sin(w_drive * t);

    % Spring Force
    k = 155; % Spring Stiffness
    F_spring = k * (x(1,1) - x_base); % Linear model
    

    % Damping Force
    F_damp = b * (x(2,1) - xd_base);

    % State evolution
    dx(1,1) = x(2,1);
    dx(2,1) = g + (F_drive - F_spring - F_damp)/m;
    dx(3,1) = xd_base;
    dx(4,1) = F_spring * (x(2,1) - xd_base);

end