%% Function to obtain the displacements of each spring in the 2-spring stack for a given force (F)

function x_eq = static_eqb(ht1_ratio, tau1, ht2_ratio, tau2, F)
    
    % Calculate the displacement of the top of spring 2
    fun2 = @(x)static_pos(x, ht2_ratio, tau2, F); 
    x2 = fsolve(fun2, 0); % Displacement of the top of spring 2

    % Calculate the displacement of the top of spring 2
    fun1 = @(x)static_pos(x, ht1_ratio, tau1, F);
    x1 = x2 + fsolve(fun1,0); % Displacement of the top of spring 1

    x_eq = [x1;x2]; % Displacements of both springs at equilibrium
    

    %{
    fun3 = @(x)static_pos(x, ht1_ratio, tau1, F);
    x3 = fsolve(fun3,0);
    x_eq = [x3;0];
    %}
end