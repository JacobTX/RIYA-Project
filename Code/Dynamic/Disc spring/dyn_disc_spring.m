function dx = dyn_disc_spring(t,x,m,g,b,ht_ratio,tau)

    % Base displacement
    x_base_max = 0;
    f_base = 10;
    w_base = 2*pi*f_base;
    x_base = x_base_max * sin(w_base * t);
    xd_base = x_base_max * w_base * cos(w_base * t);

    % Driving Force
    F0 = 0;
    w_drive = 10;
    F_drive = F0 * sin(w_drive * t);

    % Spring Force
    % k = 1
    %F_spring = k * (x(1,1) - x_base); % Linear model
    F_spring = disc_spring_force(x(1,1) - x_base,ht_ratio,tau); % Disc-spring model

    % Damping Force
    F_damping = b * (x(2,1) - xd_base);

    % State evolution
    dx(1,1) = x(2,1);
    dx(2,1) = g + (F_drive - F_spring - F_damping)/m;
    dx(3,1) = xd_base;
end