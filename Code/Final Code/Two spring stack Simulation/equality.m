%% Function that represents the equality of forces in the two springs

function F = equality(y,x_top,x_bottom, ht1_ratio, tau1, ht2_ratio, tau2)
        F = disc_spring_force(x_top - y, ht1_ratio, tau1) - disc_spring_force(y - x_bottom, ht2_ratio, tau2);
end