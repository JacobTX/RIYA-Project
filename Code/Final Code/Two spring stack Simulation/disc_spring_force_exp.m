%% Function for the spring force at a given deflection based on experiemntal data

function F = disc_spring_force_exp(x)
    x_key = [0, 0.241, 0.5, 1, 1.5, 2, 2.5]; % Key values of deflection
    F_key = [0, 41.8, 80.5, 111, 116, 123, 198]; % Key values of force corresponding to the key deflection points

    %x_key = [0, 0.241, 0.5, 1, 1.5, 1.75, 2, 2.27, 2.5];
    %F_key = [0, 37, 76, 106, 110, 112, 118, 148, 192];

    F = spline(x_key,F_key,x); % Cubic spline interpolation
end