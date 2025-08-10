clear all;
close all;
clc;
Colors = lines(6);
%================================================================================
% Parameters
%================================================================================
%============================================================================
%%      Material Properties and Spring Geometry
%=============================================================================
E = 200*10^9/1e6;  % N/mm^2
a = 34.5/2;        % outer diameter, mm
b = 22.4/2;        % inner diameter, mm
t = 0.5;           % thickness, mm
% Spring 1 h/t ratios
ht1_list =1.69;%[1.18:.01:1.41];%[1.30:.04:2.2];
% Spring 2 h/t ratios
ht2_list =1.75;%[1.19:.01:1.41];%[1.28:.04:2.2];
intercepts = [3.25:.02:3.51];

% step size [mm]
dx = .04; % Smalleset variation in delta_st

% used in algorithm for detecting snap-through evenets
closeness_tolerance = .1475;

% forward or reverse sweep (1 or -1);
direction = 1;

%=============================================================
% constants (Don't Modify)
%===============================================================
n = 2; % number of springs (keep at 2)
gamma = a/b; % alpha in the paper (ratio of outer to inner diameter)
c1 = ((gamma+1)/(gamma-1)-2/log(gamma)).*t; % M in the paper
c2 = t.^3/6*log(gamma); % N in the paper
%=============================================================
tic;
% sweep through list of h1/t and h2/t
for p = 1:1:length(ht1_list);
for o = [1:1:length(ht2_list)]
%     ht2_list = -ht1_list(p) + intercepts(o);
%     h = t.*[ht1_list(p) -ht1_list(p) + intercepts(o)];%;h = t.*[ht1_list(p) ht2_list(o)];
% x_stack = [.5:dx:(ht1_list(p)+ (-ht1_list + intercepts(o))*.9)];
    h = t.*[ht1_list(p) ht2_list(o)];
    x_stack = [.2:dx:(ht1_list(p)+ (ht2_list(o))*.95)]; % Displacement of the top of the stack (= delta_st in the paper)

%(ht1_list(p)+ (ht2_list(o))*.95) - estimate of total stroke 
        

for i = 1:1:length(x_stack);
%===========================================================================================
%=================================  VPA SOLVE ===================================================================
%===================================================================================================================
    x = sym('x',[1 n]); % Create symbolic variable x
    syms x_stack_sym; % Create symbolic variable x_stack_sym
    
    % Case for 2-spring stack
    % stack_equations_mm(1) is the expression for force in the stack
    % stack_equation_mm(2:n) represent the expressions for forces being
    % equal in consecutive disc-springs
    if n == 2;
    stack_equations_mm(1) =  E*pi/a^2.*(gamma/(gamma-1))^2.*(.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(n)*c1*(x_stack(i)-x(n-1)).^2 + (h(n)^2*c1 + c2).*(x_stack(i)-x(n-1))) -x(n);
    stack_equations_mm(n) =  (.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(2)*c1*(x_stack(i)-x(n-1)).^2 + (h(2)^2*c1 + c2).*(x_stack(i)-x(n-1)))-...
    (.5*c1.*x(1).^3 - 3/2*h(1)*c1*x(1).^2 + (h(1)^2*c1 + c2).*x(1));

% elseif n == 3;
% stack_equations_mm(1) =  E*pi/a^2.*(gamma/(gamma-1))^2.*(.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(n)*c1*(x_stack(i)-x(n-1)).^2 + (h(n)^2*c1 + c2).*(x_stack(i)-x(n-1))) -x(n);
% stack_equations_mm(2) =  (.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(3)*c1*(x_stack(i)-x(n-1)).^2 + (h(3)^2*c1 + c2).*(x_stack(i)-x(n-1)))-...
%     (.5*c1.*(x(n-1)-x(n-2)).^3 - 3/2*h(2)*c1*(x(n-1)-x(n-2)).^2 + (h(2)^2*c1 + c2).*(x(n-1)-x(n-2)));
% stack_equations_mm(n) =  (.5*c1.*(x(2)-x(1)).^3 - 3/2*h(2)*c1*(x(2)-x(1)).^2 + (h(2)^2*c1 + c2).*(x(2)-x(1)))-...
%     (.5*c1.*x(1).^3 - 3/2*h(1)*c1*x(1).^2 + (h(1)^2*c1 + c2).*x(1));
% 
% elseif n > 3;
%         stack_equations_mm(1) =  E*pi/a^2.*(gamma/(gamma-1))^2.*(.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(n)*c1*(x_stack(i)-x(n-1)).^2 + (h(n)^2*c1 + c2).*(x_stack(i)-x(n-1))) -x(n);
%         stack_equations_mm(2) =  (.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(n)*c1*(x_stack(i)-x(n-1)).^2 + (h(n)^2*c1 + c2).*(x_stack(i)-x(n-1)))-...
%     (.5*c1.*(x(n-1)-x(n-2)).^3 - 3/2*h(3)*c1*(x(n-1)-x(n-2)).^2 + (h(3)^2*c1 + c2).*(x(n-1)-x(n-2)));
%         
%     for i = [3:1:n-1];
%     
%     stack_equations(i) =  (.5*c1.*(x(i)-x(i-1)).^3 - 3/2*h(i)*c1*(x(i)-x(i-1)).^2 + (h(i)^2*c1 + c2).*(x(i)-x(i-1)))-...
%     (.5*c1.*(x(i-1)-x(i-2)).^3 - 3/2*h(i-1)*c1*(x(i-1)-x(i-2)).^2 + (h(i-1)^2*c1 + c2).*x(i-1)-x(i-2));
%     
%     end
% stack_equations_mm(n) =  (.5*c1.*(x(2)-x(1)).^3 - 3/2*h(2)*c1*(x(2)-x(1)).^2 + (h(2)^2*c1 + c2).*(x(2)-x(1)))-...
%     (.5*c1.*x(1).^3 - 3/2*h(1)*c1*x(1).^2 + (h(1)^2*c1 + c2).*x(1));
% 
    end
%     

% Solve for stack force and displacement x1 (= delta_1 in the paper)
 for u = 1:n;
     assume(x(u),'real');
 end;
    X = vpasolve(stack_equations_mm);  % Solve simultaneous nonlinear equations
    fname = ['x',num2str(n)];
    forceVPA = X.(fname);                % solution for stack force
    X.x1;                               % solution for x1
    x_st_name = ['x_st',num2str(i)];
    X_full_sol.(x_st_name) = X;
    
    
%
%=======================================================================================================
%%%======================================= Data Processing and Analysis%%%================================
%===========================================================================================================

% Plot each plot in a different gray color.
M = zeros(n+1,1);
x1_varied = [.2:.02:2.1];  % varying x1 to calulate potential energy at a given step
ramp = linspace(0, 0.75, length(x_stack));
listOfGrayColors = [ramp; ramp; ramp]';

    for j = 1:length(forceVPA);  % sweep through all solutions at a given step of x_stack 
        for k = 1;  % set k=1 if you want to plot force-deflection curves, set k = 1:1:length(x1_varied) if you want to plot potential energy maps
        if n == 2;   % n = 2 always (n= number of springs in stack)
            syms x1;
%         force_vector = .6*min(forceVPA):.1:1.4*max(forceVPA);
%         V1(j) = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*X.x1(j)^4 - 1/2*h(1)*c1*X.x1(j)^3 + (h(1)^2*c1+c2)/2*X.x1(j)^2);  % spring 1 potential energy
%         V2(j) = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*(x_stack(i)-X.x1(j))^4 - 1/2*h(2)*c1*(x_stack(i)-X.x1(j))^3 + (h(2)^2*c1+c2)/2*(x_stack(i)-X.x1(j))^2);  % spring 2 potential energy
%         Vtot1(j) = V1(j) + V2(j);
%         Vtot(j) = V1(j) + V2(j)- forceVPA(j).*x_stack(i);
%         Vtot_curve = V1(j) + V2(j) - force_vector.*x_stack(i);
%===========================================================================================================  
%%%%========================== Symbolic potential energy expressions =====================================
%===================================================================================================================
        V1_sym = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*x1^4 - 1/2*h(1)*c1*x1^3 + (h(1)^2*c1+c2)/2*x1^2);  % spring 1 potential energy
        V2_sym = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*(x_stack_sym-x1)^4 - 1/2*h(2)*c1*(x_stack_sym-x1)^3 + (h(2)^2*c1+c2)/2*(x_stack_sym-x1)^2);  % spring 2 potential energy
        Vtot_sym = V1_sym + V2_sym;%- forceVPA(j).*x_stack_sym;
        Veval(j) = subs(Vtot_sym,[x1 x_stack_sym],[X.x1(j) x_stack(i)]);
        Veval(i) = subs(Vtot_sym,[x1],[1.0803]);
        % Plot the potential energy curve
        %x_stack_vary = .89:.005:3.5;
        %for z = 1:length(x_stack_vary)
             %Veval2(z) = subs(Vtot_sym,[x1 x_stack_sym],[X.x1(j) x_stack_vary(z)]); % Evaluate the value of total system potential energy
        %end
        %figure(10);
        %plot(x_stack_vary,Veval2); % Obtain the potential energy curve as a function of stack displacement
         
%%%===========================================================================================
%%======================== For Surface Plot Generation===================================
%%%============================================================================================
        % Potential Energy while varying x1
         V1_sym = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*x1_varied(k)^4 - 1/2*h(1)*c1*x1_varied(k)^3 + (h(1)^2*c1+c2)/2*x1_varied(k)^2);  % spring 1 potential energy
         V2_sym = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*(x_stack_sym-x1_varied(k))^4 - 1/2*h(2)*c1*(x_stack_sym-x1_varied(k))^3 + (h(2)^2*c1+c2)/2*(x_stack_sym-x1_varied(k))^2);  % spring 2 potential energy
         calculated_force1 = E*pi/a^2.*(gamma/(gamma-1))^2.*(.5*c1.*(x_stack_sym-x1_varied(k)).^3 - 3/2*h(n)*c1*(x_stack_sym-x1_varied(k)).^2 + (h(n)^2*c1 + c2).*(x_stack_sym-x1_varied(k)));
         calculated_force2 = E*pi/a^2.*(gamma/(gamma-1))^2.*(.5*c1.*x1_varied(k).^3 - 3/2*h(1)*c1*x1_varied(k).^2 + (h(1)^2*c1 + c2).*x1_varied(k));
         CF1(k) = eval(subs(calculated_force1,[x_stack_sym],[x_stack(i)]));
         CF2(k) = eval(subs(calculated_force2,[x_stack_sym],[x_stack(i)]));
         
         Vtot_sym_surf = V1_sym + V2_sym;%- CF1(k).*x_stack_sym;
         Veval(k) = subs(Vtot_sym_surf,[x_stack_sym],[x_stack(i)]);
%%%==========================================================================================
%%==========================Determine Stability of Each Point ===================================
%%%============================================================================================
               
        % first derivative with respect to all variables
              
        dP_dx1 = diff(Vtot_sym,x1); % Evaluate the first derivative of potential

%         dP_dxst = diff(Vtot_sym,x_stack_sym);
        
%         test = subs(dP_dx1,[x1 x_stack_sym],[X.x1(j) x_stack(i)])
        
        d2P_dx12 = diff(dP_dx1,x1); % Expression for second derivative of potential
%         d2P_dx1dxst = diff(dP_dx1,x_stack_sym);
%         d2P_dxstdx1 = diff(dP_dxst,x1);
%         d2P_dxst2 = diff(dP_dxst,x_stack_sym);
%         
%         
%         H = [d2P_dx12 d2P_dx1dxst; d2P_dxstdx1 d2P_dxst2];
%         H_eval = subs(H,[x1 x_stack_sym],[X.x1(j) x_stack(i)]);
% %         Ht = [H_eval(1,1) H_eval(1,2); H_eval(2,1) H_eval(2,2)];
%        EV = eig(H_eval);
       d2P_dx12_eval = subs(d2P_dx12,[x1 x_stack_sym],[X.x1(j) x_stack(i)]); % Evaluate the second derivative of the potential for different values of x1 and x_st
        
%==================================================================
%Plotting each solution
%===================================================================
        if d2P_dx12_eval > 0;  %all(EV(:)>0) == 1; % Stable equilibrium
       figure(o);
%        yyaxis left
       p1 = plot(x_stack(i),forceVPA(j),'o','Color',Colors(1,:)); 
       p1.MarkerFaceColor = Colors(1,:);
       hold on;
       grid on;
       set(gca,'ylim',[60 160],'FontSize',20);
%        yyaxis right
%        p2 = plot(x_stack(i),X.x1(j),'^','Color',Colors(1,:));
%        p2.MarkerFaceColor = Colors(1,:);
%        hold on;
       
        stability{i,j} = 1;
        elseif d2P_dx12_eval < 0; % Unstable equilibrium
        figure(o);
%         yyaxis left
        p1 = plot(x_stack(i),forceVPA(j),'o','Color',Colors(1,:));
        hold on;
        grid on;
%         yyaxis right
%         plot(x_stack(i),X.x1(j),'^','Color',Colors(1,:));
%         hold on;
         stability{i,j} = 0;
        else
%             fprintf('error');
        end
%         
% %           M(1:n+1,j) = [forceVPA(j); EV];
% %         
% %         elseif n ==3;
% %         else n > 3;
% %         V2(j) = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*(x2(j)-x1(j))^4 - 1/2*h(2)*c1*(x2(j)-x1(j))^3 + (h(2)*c1+c2)/2*(x2(j)-x1(j))^2);  % spring 2 potential energy
% %         V3(j) = E*pi/a^2*(gamma/(gamma-1))^2.*(1/8*c1*(xstack(i)-x2(j))^4 - 1/2*h(3)*c1*(xstack(i)-x2(j))^3 + (h(3)*c1+c2)/2*(xstack(i)-x2(j))^2);  % spring 3 potential energy
% %         Vtot(j) = V1(j) + V2(j) + V3(j) - forceVPA(j).*x_stack(i);
        end
        
        clear EV;
        
        %%%==================================================================================================================
    end
%      [Vmin Vminidx] = min(Vtot1);
%     force_min(i) = forceVPA(Vminidx);
%     figure(1);
%     p1 = plot(x_stack(i).*ones(1,length(forceVPA)),forceVPA,'o','Color',lines(1));
% % % p1 = plot(x_stack(i),force_min(i),'o','Color',lines(1));
%     p1.MarkerSize = 5;
%     hold on;
%     plot(x_stack(i),CF1,'o',x_stack(i),CF2,'^');
%     hold on;
%     [forceVPA Veval'];
%     M = zeros(n+1,1);
%     clear force_min;
               

    end
    

%     Veval_surf(i,1:length(x1_varied)) = Veval./max(Veval);
%     
% figure(4);
% title(strcat('x_stack = ',num2str(x_stack(i))));
% plot(x1_varied, Veval,'Color', listOfGrayColors(i,:),'LineWidth',1.5);
% hold on;
% set(gca,'FontSize',18,'xlim',[.3 1.4],'ylim',[80 240],'Position',[   0.1939    0.1467    0.7111    0.7632]);
% xlabel('\delta_1 [mm]');
% ylabel('U [N-mm]')

% figure(10+i);
% plot(x1_varied,CF1,x1_varied,CF2);
    clear Vtot1;
   clear H_eval
   clear H
    clear Vtot;
    clear Veval;
    clear Vtot_curve;
    clear CF1;
    clear CF2;
    %======================================================================================================
%     clear forceVPA; clear x1; clear x2; clear Vtot
end

%      figure(4);
% surf(x1_varied,x_stack,double(Veval_surf));
%  xlabel('x1 [mm]');
% ylabel('Stack Displacement [mm]')
% zlabel('Potential Energy');
% hold on; 
toc
%%
% Code for plotting forward and reverse force-deflection curves and
% calculation stack stiffness from those curves
clc
clear X;
syms d;

if direction ==1;
    steps = 1:1:length(x_stack);
elseif direction == -1;
    steps = length(x_stack):-1:1;
end
for i = steps
    % current solution
    x_st_name = ['x_st',num2str(i)];
    X = X_full_sol.(x_st_name);
    forceVPA = X.(fname);
    x1_solution = X.x1;
     
    idpath = 1;
    if length(x1_solution)==1;
        x1path(i) = x1_solution;
        fpath(i) = forceVPA;
       
    else  % If more than 1 solution exists:
        
         %  same number of equilibrium point?
        if length(x1_solution) == length(x1_prev_all);
            num_points_created = 0;
            num_points_destroyed = 0;
        elseif length(x1_solution) < length(x1_prev_all);  % points destroyed
            num_points_created = 0;
            num_points_destroyed = length(x1_prev_all)-length(x1_solution); 
        elseif length(x1_solution) > length(x1_prev_all); %points created
            num_points_created = length(x1_solution)-length(x1_prev_all);
            num_points_destroyed = 0;
        end
        
     stable_solutions_idx = find(cell2mat(stability(i,:))==1);
     %=========================================================================================================================
     %==========  Is there a new stable fixed point near the previous point and no points destroyed? %===========================
     %===========================================================================================================================
%      for z = 1:length(stable_solutions_idx);
         if direction == 1;
             % distance between previous point and new stable points
         closeness_to_previous_pathpoint = abs(x1path(i-1)-x1_solution(stable_solutions_idx));
         [min_closeness min_closeness_idx] = min(closeness_to_previous_pathpoint);
         
         % distance between previous point and previous unstable points
         previous_unstable_solutions_idx = find(cell2mat(stability(i-1,:))==0);
         previous_point_closeness_to_all_previous_unstable_points = abs(x1path(i-1)-x1_prev_all(previous_unstable_solutions_idx));
         else direction == -1;
              % distance between previous point and new stable points
             closeness_to_previous_pathpoint = abs(x1path(i+1)-x1_solution(stable_solutions_idx));
            [min_closeness min_closeness_idx] = min(closeness_to_previous_pathpoint);
            
            % distance between previous point and previous unstable points
            previous_unstable_solutions_idx = find(cell2mat(stability(i+1,:))==0);
            previous_point_closeness_to_all_previous_unstable_points = abs(x1path(i+1)-x1_prev_all(previous_unstable_solutions_idx));
         end;
         
          
          
         if min_closeness <= closeness_tolerance & num_points_destroyed == 0;
             x1path(i) = x1_solution(stable_solutions_idx(min_closeness_idx));
             fpath(i) = forceVPA(stable_solutions_idx(min_closeness_idx));
             
         elseif min_closeness >= closeness_tolerance & num_points_destroyed == 0;
             fprintf('No stable fixed point found near previous point and no points destroyed \n');
             [i h]
             
         elseif min_closeness >=closeness_tolerance & num_points_destroyed > 0;
              % find gradient at estimated bifurcation point and current x_stack; then method of steepest descnet to find next stable x_path
              % evaluate gradient at first point:  this is the direction
             fprintf('Min closeness greater than 2dx and points destroyed');
         end 
    end
    
    %Determining if a snap-through event has occured at current step
    if i == 1 | i==length(x_stack);
        snapthrough_set(i) = 0;  % no snap-through at initial point (x=0 or x=x_final)
    else
        x1_diff(i) = abs(x1path(i)-x1_path_prev);
        if abs(x1path(i)-x1_path_prev) >= closeness_tolerance;  % snapthrough occurs if there is a large jump in x1 from previous step
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
 
 % Determining if a snap-through event has occured over the whole
 % force-deflection curve for a given h1/t and h2/t
 if sum_snap_through == 0;
     snapthrough_matrix(o,p) = 0;
     Colors = [.85 .85 .85];
 elseif sum_snap_through >= 1;
     snapthrough_matrix(o,p) =1;
     Colors =  [0 0 0];
 end

 
 clear x1_diff;
% figure(o);
% hold on;
%     yyaxis left
%     plot(x_stack,fpath,'-','Color',[0 0 0],'LineWidth',1.5);
%     hold on;
%     set(gca,'FontSize',18);
%     grid on;
%     yyaxis right;
%     plot(x_stack,x1path,'-','Color',[0 0 0],'LineWidth',1.5);
%     set(gca,'FontSize',18);
%     grid on;

%Calculating stack stiffness
    for k = 1:length(fpath)-1;
        stiffness(k) = (fpath(k+1)-fpath(k))/dx;
    end
    min_stiffness(o,p) = min(stiffness);
%     plot(x_stack(1:end-1),stiffness,'o-','LineWidth',1.4);
%     set(gca,'FontSize',18);
%     grid on;
    clear stiffness;
% figure(13);
% plot(ht1_list(p),ht2_list(o),'s','Color',Colors,'MarkerSize',14,'MarkerFaceColor',Colors);
% hold on;
clear x1path;
clear fpath;
clear stability;
end
end
toc