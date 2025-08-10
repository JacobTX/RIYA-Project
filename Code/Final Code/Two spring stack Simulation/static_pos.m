%% Function used to solve for the spring displacement (x) for a given force (F)

function S = static_pos(x, ht_ratio, tau, F)
       S = disc_spring_force(x, ht_ratio, tau) - F;
       %S = disc_spring_force_exp(x)-F;
end