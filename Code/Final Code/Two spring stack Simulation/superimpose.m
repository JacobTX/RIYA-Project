%% Physical domain behaviour

%% Plot the force vs displacement behaviour in dynamic situation
fig = figure;
set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
plot(x_rel, fs1, 'LineWidth', 2, 'Color', 'Blue')
set(gca,'FontSize',15);
hold on


%% Obtain static force-deflection curve for 2-spring stack and superimpose

% The code below is used to solve the 2-spring stack non-linear equations
% in order to obtain the force-deflection behaviour and potential energy
% curves of the 2-spring stack

Colors = lines(6);

%% Parameters - Material Properties and Spring Geometry

E = 200*10^9/1e6;  % N/mm^2
a = 34.5/2;        % outer diameter, mm
b = 22.4/2;        % inner diameter, mm
t = 0.5;           % thickness, mm


ht1_list = ht1_ratio; % h/t ratio of Spring 1
ht2_list = ht2_ratio; % h/t ratio of Spring 2

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
    x_stack = .1:dx:x_max; % Displacement of the top of the stack (= delta_st in the paper)


%% Iterate over increasing values of x_stack

    for i = 1:1:length(x_stack)
        x = sym('x',[1 n]); % Create symbolic vector variable x with dimensions 1 x n
        syms x_stack_sym; % Create symbolic variable x_stack_sym
        
        % Case for 2-spring stack
        % stack_equations_mm(1) is the expression for force in the stack
        % stack_equation_mm(2:n) represent the expressions for forces being
        % equal in consecutive disc-springs

        stack_equations_mm(1) =  E*pi/a^2.*(gamma/(gamma-1))^2.*(.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(n)*c1*(x_stack(i)-x(n-1)).^2 + (h(n)^2*c1 + c2).*(x_stack(i)-x(n-1))) -x(n);
        stack_equations_mm(2) =  (.5*c1.*(x_stack(i)-x(n-1)).^3 - 3/2*h(2)*c1*(x_stack(i)-x(n-1)).^2 + (h(2)^2*c1 + c2).*(x_stack(i)-x(n-1)))-(.5*c1.*x(1).^3 - 3/2*h(1)*c1*x(1).^2 + (h(1)^2*c1 + c2).*x(1));

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
                       p1 = plot(x_stack(i),forceVPA(j),'o','Color',Colors(1,:),'DisplayName','Static');
                       p1.MarkerFaceColor = Colors(1,:);
                       hold on;
                       grid on;
                       set(gca,'ylim',[40 250],'FontSize',15);
                       title('Force-deflection curve','Interpreter','latex','FontSize',15)
                       xlabel('$x_{rel}$ [mm]','Interpreter','latex','FontSize',15)
                       ylabel('Force [N]','Interpreter','latex','FontSize',15)
                       set(gca,'FontSize',15);
       

                    elseif d2P_dx12_eval < 0 % Unstable equilibrium
                       %% === Force vs x_stack ===
                       set(gcf,'defaultAxesTickLabelInterpreter','latex'); 
                       p1 = plot(x_stack(i),forceVPA(j),'+','Color',Colors(1,:));
                       set(gca,'FontSize',15);
                       hold on;

                    else
                      fprintf('error');

                    end       
         end
     end
end 
end

hold on
scatter(x_rel(1),fs1(1),100,'red','filled')
legend("Dynamic","Static",'Interpreter','latex','FontSize',15)

%exportgraphics(fig,'physdom.pdf','Resolution',300)