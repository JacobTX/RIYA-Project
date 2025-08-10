%% Annotated and simplified version of the MATLAB Code used to analyse the 2-spring stack

% The code below is used to solve the 2-spring stack non-linear equations
% in order to obtain the force-deflection behaviour and potential energy
% curves of the 2-spring stack

clear all;
close all;
clc;
Colors = lines(6);
%% Parameters - Material Properties and Spring Geometry

E = 200*10^9/1e6;  % N/mm^2
a = 34.5/2;        % outer diameter, mm
b = 22.4/2;        % inner diameter, mm
t = 0.5;           % thickness, mm

% h/t ratios of Spring 1
ht1_list = 1.4;
%ht1_list =1.69;
%ht1_list = 1.18:.01:1.41;
%ht1_list = 1.30:.04:2.2;

% h/t ratios of Spring 2
ht2_list = 1.4;
%ht2_list =1.75;
%ht2_list = 1.19:.01:1.41;
%ht2_list = 1.28:.04:2.2;

dx = .04; % step size [mm] - Smalleset variation in delta_st
closeness_tolerance = .1475; % Used in algorithm for detecting snap-through evenets
direction = 1; % Forward or reverse sweep (1 or -1);

%% Constants (Don't Modify)

n = 2; % Number of Springs (keep at 2)
gamma = a/b; % alpha in the paper (ratio of outer to inner diameter)
c1 = ((gamma+1)/(gamma-1)-2/log(gamma)).*t; % M in the paper
c2 = t.^3/6*log(gamma); % N in the paper

%% Sweep through list of h1/t and h2/t to analyze for each pair of (h1/t, h2/t)

for p = 1:1:length(ht1_list)
for o = 1:1:length(ht2_list)

    h = t.*[ht1_list(p) ht2_list(o)]; % Heights [h1,h2] of the springs 1 and 2 in the stack
    x_max = (ht1_list(p)+ (ht2_list(o))*.95); % Conservative estimate of total stroke 
    x_stack = .2:dx:x_max; % Displacement of the top of the stack (= delta_st in the paper)


%% Iterate over increasing values of x_stack

    for i = 1:1:length(x_stack)
        x = sym('x',[1 n]); % Create symbolic vector variable x with dimensions 1 x n
        syms x_stack_sym; % Create symbolic variable x_stack_sym
        
        % Case for 2-spring stack
        % stack_equations_mm(1) is the expression for force in the stack
        % stack_equation_mm(2:n) represent the expressions for forces being
        % equal in consecutive disc-springs
    
        if n == 2
            stack_equations_mm(1) =  E*pi/a^2.*(gamma/(gamma-1))^2.*(.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(n)*c1*(x_stack(i)-x(n-1)).^2 + (h(n)^2*c1 + c2).*(x_stack(i)-x(n-1))) -x(n);
            stack_equations_mm(n) =  (.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(2)*c1*(x_stack(i)-x(n-1)).^2 + (h(2)^2*c1 + c2).*(x_stack(i)-x(n-1)))-(.5*c1.*x(1).^3 - 3/2*h(1)*c1*x(1).^2 + (h(1)^2*c1 + c2).*x(1));
        end
        
        % Solve for stack force and displacement x1 (= delta_1 in the paper)
    
        for u = 1:n % Sweep through each displacement component 
            assume(x(u),'real'); % Letting x(u) be a real value only
        end

        % Solve simultaneous nonlinear equations using vpasolve
        X = vpasolve(stack_equations_mm);  
        fname = ['x',num2str(n)]; % fname = 'x2'
        forceVPA = X.(fname);   % For n = 2, forceVPA = X.x2             
        x_st_name = ['x_st',num2str(i)];
        X_full_sol.(x_st_name) = X;
    
        %% Data Processing and Analysis
       
        for j = 1:length(forceVPA)  % Sweep through all solutions at a given step of x_stack 
            % Set k=1 if you want to plot force-deflection curves
            for k = 1
                if n == 2   % n = 2 in this case (n = number of springs in stack)
                    syms x1;

                    % Symbolic potential energy expressions
                    V1_sym = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*x1^4 - 1/2*h(1)*c1*x1^3 + (h(1)^2*c1+c2)/2*x1^2);  % Spring 1 Potential Energy
                    V2_sym = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*(x_stack_sym-x1)^4 - 1/2*h(2)*c1*(x_stack_sym-x1)^3 + (h(2)^2*c1+c2)/2*(x_stack_sym-x1)^2);  % Spring 2 potential energy
                    Vtot_sym = V1_sym + V2_sym; % Total potential energy of the stack = - forceVPA(j).*x_stack_sym;
                  
                    %% Determine Stability of Each Point       
                    dP_dx1 = diff(Vtot_sym,x1); % Evaluate the first derivative of potential
                    d2P_dx12 = diff(dP_dx1,x1); % Expression for second derivative of potential
                    d2P_dx12_eval = subs(d2P_dx12,[x1 x_stack_sym],[X.x1(j) x_stack(i)]); % Evaluate the second derivative of the potential for different values of x1 and x_st
        

                    % Plotting each solution

                    if d2P_dx12_eval > 0  % Stable equilibrium
                       %% === Force vs x_stack ===
                       
                       set(gcf,'position',[0,0,1500,750]);
                       subplot(1,3,1);
                       p1 = plot(x_stack(i),forceVPA(j),'o','Color',Colors(1,:)); 
                       p1.MarkerFaceColor = Colors(1,:);
                       hold on;
                       grid on;
                       set(gca,'ylim',[0 200],'FontSize',15);
                       title('Force-deflection curve')
                       xlabel('x_{stack} [mm]')
                       ylabel('Force [N]')
                       set(gca,'FontSize',15);

                       %% === Potential Energy vs x_stack ===
                       Veval = subs(Vtot_sym,[x1 x_stack_sym],[X.x1(j) x_stack(i)]); % Evaluate the value of total system potential energy
                    
                       subplot(1,3,2);
                       plot(x_stack(i),Veval,'o','Color','Red');
                       hold on
                       title('Potential Energy vs x_{stack}')
                       xlabel('x_{stack} [mm]')
                       ylabel('Potential Energy')
                       set(gca,'FontSize',15);
                       
                       %% === x1 vs x_stack ===
                       subplot(1,3,3);
                       p2 = plot(x_stack(i),X.x1(j),'^','Color',Colors(1,:));
                       p2.MarkerFaceColor = Colors(1,:);
                       hold on;
                       grid on;
                       set(gca,'FontSize',15);
                       title('x_1 vs x_{stack}')
                       xlabel('x_{stack} [mm]')
                       ylabel('x_1 [mm]')
                       set(gca,'FontSize',15);
       
                       stability{i,j} = 1;

                    elseif d2P_dx12_eval < 0 % Unstable equilibrium
                       %% === Force vs x_stack ===
                       subplot(1,3,1);
                       p1 = plot(x_stack(i),forceVPA(j),'+','Color',Colors(1,:));
                       hold on;

                       %% === Potential Energy vs x_stack ===
                       Veval = subs(Vtot_sym,[x1 x_stack_sym],[X.x1(j) x_stack(i)]); % Evaluate the value of total system potential energy
                       subplot(1,3,2);
                       plot(x_stack(i),Veval,'+','Color','Black'); % Obtain the potential energy curve as a function of stack displacement
                       hold on

                       %% === x1 vs x_stack ===
                       subplot(1,3,3);
                       p2 = plot(x_stack(i),X.x1(j),'^','Color',Colors(1,:));
                       hold on;
                  
                       stability{i,j} = 0;

                    else
                      fprintf('error');

                    end       

                end
            end
        end
    end 

%% Code for plotting forward and reverse force-deflection curves and calculation stack stiffness from those curves
clc
clear X;
syms d;

if direction ==1
    steps = 1:1:length(x_stack);
elseif direction == -1
    steps = length(x_stack):-1:1;
end

for i = steps
    % current solution
    x_st_name = ['x_st',num2str(i)];
    X = X_full_sol.(x_st_name);
    forceVPA = X.(fname);
    x1_solution = X.x1;
     
    idpath = 1;
    if length(x1_solution)==1
        x1path(i) = x1_solution;
        fpath(i) = forceVPA;
       
    else  % If more than 1 solution exists:
        
         %  same number of equilibrium point?
        if length(x1_solution) == length(x1_prev_all)
            num_points_created = 0;
            num_points_destroyed = 0;
        elseif length(x1_solution) < length(x1_prev_all)  % points destroyed
            num_points_created = 0;
            num_points_destroyed = length(x1_prev_all)-length(x1_solution); 
        elseif length(x1_solution) > length(x1_prev_all) %points created
            num_points_created = length(x1_solution)-length(x1_prev_all);
            num_points_destroyed = 0;
        end
        
     stable_solutions_idx = find(cell2mat(stability(i,:))==1);
     %=========================================================================================================================
     %==========  Is there a new stable fixed point near the previous point and no points destroyed? %===========================
     %===========================================================================================================================
%      for z = 1:length(stable_solutions_idx);
         if direction == 1
             % distance between previous point and new stable points
            closeness_to_previous_pathpoint = abs(x1path(i-1)-x1_solution(stable_solutions_idx));
            [min_closeness, min_closeness_idx] = min(closeness_to_previous_pathpoint);
         
            % distance between previous point and previous unstable points
            previous_unstable_solutions_idx = find(cell2mat(stability(i-1,:))==0);
            previous_point_closeness_to_all_previous_unstable_points = abs(x1path(i-1)-x1_prev_all(previous_unstable_solutions_idx));

         elseif direction == -1
            % distance between previous point and new stable points
            closeness_to_previous_pathpoint = abs(x1path(i+1)-x1_solution(stable_solutions_idx));
            [min_closeness, min_closeness_idx] = min(closeness_to_previous_pathpoint);
            
            % distance between previous point and previous unstable points
            previous_unstable_solutions_idx = find(cell2mat(stability(i+1,:))==0);
            previous_point_closeness_to_all_previous_unstable_points = abs(x1path(i+1)-x1_prev_all(previous_unstable_solutions_idx));
         end
        
         if min_closeness <= closeness_tolerance && num_points_destroyed == 0
             x1path(i) = x1_solution(stable_solutions_idx(min_closeness_idx));
             fpath(i) = forceVPA(stable_solutions_idx(min_closeness_idx));
             
         elseif min_closeness >= closeness_tolerance && num_points_destroyed == 0
             fprintf('No stable fixed point found near previous point and no points destroyed \n');
             
         elseif min_closeness >=closeness_tolerance && num_points_destroyed > 0
              % find gradient at estimated bifurcation point and current x_stack; then method of steepest descnet to find next stable x_path
              % evaluate gradient at first point:  this is the direction
             fprintf('Min closeness greater than 2dx and points destroyed');
         end 
    end
 
    %Determining if a snap-through event has occured at current step
    if i == 1 || i==length(x_stack)
        snapthrough_set(i) = 0;  % no snap-through at initial point (x=0 or x=x_final)
    else
        x1_diff(i) = abs(x1path(i)-x1_path_prev);
        if abs(x1path(i)-x1_path_prev) >= closeness_tolerance % snapthrough occurs if there is a large jump in x1 from previous step
            snapthrough_set(i) = 1;
            force_snapthrough(o,p) = fpath(i)-fprev;
        else
            snapthrough_set(i) = 0;
        end
    end
        
    x1_prev_all = x1_solution;
    x1_path_prev = x1path(i);
    fprev = fpath(i);
        
end


sum_snap_through = sum(snapthrough_set);
clear snapthrough_set;
 
%% Determining if a snap-through event has occured over the whole

if sum_snap_through == 0
     snapthrough_matrix(o,p) = 0;
     Colors = [.85 .85 .85];
elseif sum_snap_through >= 1
     snapthrough_matrix(o,p) =1;
     Colors =  [0 0 0];
end

 
clear x1_diff;

%% === Plots of Force vs x_stack and x1 vs x_stack for a given direction === 
% figure(o);
% hold on;
% yyaxis left
% plot(x_stack,fpath,'-','Color',[0 0 0],'LineWidth',1.5);
% hold on;
% set(gca,'FontSize',18);
% grid on;
% yyaxis right;
% plot(x_stack,x1path,'-','Color',[0 0 0],'LineWidth',1.5);
% set(gca,'FontSize',18);
% grid on;

%% Calculating stack stiffness
 for k = 1:length(fpath)-1
     stiffness(k) = (fpath(k+1)-fpath(k))/dx;
 end
    
 min_stiffness(o,p) = min(stiffness);

 %% === Plot of Stiffness vs x_stack ===

 %plot(x_stack(1:end-1),stiffness,'o-','LineWidth',1.4);
 %xlabel("x_{stack}")
 %ylabel("Stiffness")
 %title("Stiffness vs x_{stack}")
 %set(gca,'FontSize',15);
 %grid on;

clear stiffness;
clear x1path;
clear fpath;


end
end