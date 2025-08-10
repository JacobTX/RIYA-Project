%% Read Simulation and Experimental Acceleration Data for Experiment 1

for i=9:9 % Range from 1 to 12

    close all;

    f = 4 + 2 * (i - 1); % Frequency of base excitation
    
    % Read data from csv files in appropriate folders (change folder path accordingly)
    T1 = readtable('/Users/jacobsony/Desktop/try1/a_sim_' + string(f) + '.csv');
    T2 = readtable('/Users/jacobsony/Desktop/try2/a_sim_' + string(f) + '.csv');
    T3 = readtable('/Users/jacobsony/Desktop/try3/a_sim_' + string(f) + '.csv');
    T4 = readtable('/Users/jacobsony/Documents/IITB/RIYA/RIYA Project/Code/Dynamic/Final/Exp1 Data Analysis/set' + string(i) + '.xlsx','Sheet','Vibration');
    
    % Acceleration data for c = 0.5, mu = 0
    t1 = T1.Var1;
    a1_m = T1.Var6;
    a1_b = T1.Var7;
    
    % Acceleration data for c = 1e-4, mu = 0.04
    t2 = T2.Var1;
    a2_m = T2.Var6;
    a2_b = T2.Var7;
    
    % Acceleration data for c = 1e-4, mu = 0.04
    t3 = T3.Var1;
    a3_m = T3.Var6;
    a3_b = T3.Var7;
    
    % Experimental accelerometer data
    a4_m = 9.81 * 1000 * T4.pcb5 ;
    a4_b = 9.81 * 1000 * T4.pcb4;
    t4 = 0:(1/2000):(length(a4_b)-1)/2000;
    
    t = t2;

    % Plot acceleration profiles

    num_cyc = 6;
    start_t = max(t) - (num_cyc/f);
    start_t4 = max(t4) - (num_cyc/f);
    end_t = max(t);
    min_y = min([min(a1_m),min(a1_b),min(a2_m),min(a2_b),min(a3_m),min(a3_b),min(a4_b),min(a4_m)]);
    max_y = max([max(a1_m),max(a1_b),max(a2_m),max(a2_b),max(a3_m),max(a3_b),max(a4_b),max(a4_m)]);

    F = figure;
    set(gcf,'defaultAxesTickLabelInterpreter','latex'); 

    plot(t - start_t, a1_b,'LineWidth',2,'Color','Red')
    hold on
    plot(t4 - start_t4 + 0.55*(1/f), a4_b,'LineWidth',2,'Color','#FFA500','LineStyle',':')
    hold on
    plot(t - start_t,a1_m,'LineWidth',2,'Color','Blue')
    hold on
    plot(t - start_t, a2_m,'LineWidth',2,'Color','#50C878','LineStyle',':')
    hold on
    plot(t - start_t, a3_m,'LineWidth',2,'Color','Blue','LineStyle','-.')
    hold on
    plot(t4 - start_t4 + 0.55*(1/f), a4_m,'LineWidth',2,'Color','Magenta')
    
    title("Acceleration profile","Interpreter","latex","FontSize",15)
    legend("$\ddot{x}_{base}$(Simulation)","$\ddot{x}_{base}$(Experimental)","$\ddot{x}_{st}(c = 0.5, \mu = 0)$","$\ddot{x}_{st}(c = 10^{-4}, \mu = 0.04)$","$\ddot{x}_{st}(c = 0.2, \mu = 0.03)$","$\ddot{x}_{st}$(Experimental)","Interpreter","latex","FontSize",15)
    xlim([0, end_t - start_t])
    ylim([min_y-0.1*abs(min_y),max_y+ 1.2*abs(max_y)])
    xlabel('$t$ [s]',"Interpreter","latex","FontSize",15)
    ylabel('$\ddot{x}$ [mm/s\textsuperscript{2}]',"Interpreter","latex","FontSize",15)
    set(gca, 'FontSize',15)

    % Export plot
    plotpath = fullfile('/Users/jacobsony/Desktop/pic1.pdf'); % Use appropriate folder path
    exportgraphics(F,plotpath,'Resolution',1000)

end