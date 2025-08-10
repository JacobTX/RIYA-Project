function dx = dyn_spring_stack(t,x,m1,m2,g,b1,b2)

    % Base displacement
    x_base_max = 0;
    w_base = 10;
    x_base = x_base_max * sin(w_base * t);
    xd_base = x_base_max * w_base * cos(w_base * t);

    % Driving Force
    F0 = 0;
    w_drive = 10;
    F_drive = F0 * sin(w_drive * t);

    % Spring Forces
    k1 = 155; % Spring Stiffness
    F1_spring = k1 * (x(1,1) - x(3,1)); % Linear model

    k2 = 155;
    F2_spring = k2 * (x(3,1) - x_base);
    

    % Damping Force
    F1_damping = b1 * (x(2,1) - x(4,1));
    F2_damping = b2 * (x(4,1) - xd_base);

    % State evolution
    dx(1,1) = x(2,1);
    dx(2,1) = g + (F_drive - F1_spring - F1_damping)/m1;
    dx(3,1) = x(4,1);
    dx(4,1) = (F1_spring - F2_spring + F1_damping - F2_damping)/m2;
    dx(5,1) = xd_base;
end