%% Main code 
% AE 199 Project 2
% Daniel Huang & Emily Pippin
clc;close all;clear
%% Setup 1: y' = x*y in-class example
% Tunable parameters
h1 = .01;
x0_1 = 0;
xf_1 = 3;
y0_1 = 1;
% EoM
eqn1 = @(x,y)(x*y);
% Numerical Solution
tspan_1 = [x0_1, xf_1];
sol_1RK = RK4(eqn1,tspan_1,y0_1,h1);
sol_1ABM = ABM4(eqn1,tspan_1,y0_1,h1);
% Extract solutions
x1RK = sol_1RK.x;
y1RK = sol_1RK.y;
x1ABM = sol_1ABM.x;
y1ABM = sol_1ABM.y;
% Exact solution
y1_exactRK = exp(x1RK.^2/2);
y1_exactABM = exp(x1ABM.^2/2);
% Error
er1RK = -y1_exactRK + y1RK;
er1ABM = -y1_exactABM + y1ABM;
% Plotting
figure
subplot(2,1,1)
plot(x1RK,y1RK,x1ABM,y1ABM,'-.',x1RK,y1_exactRK,'--','linewidth',3)
title('RK4 vs ABM4 vs Exact solution')
xlabel('x');ylabel('y')
legend('RK4','ABM4','Exact solution','Location','best')
subplot(2,1,2)
plot(x1RK,er1RK,x1ABM,er1ABM,'linewidth',2.5)
hold on
plot(tspan_1,[0 0],'--k')
hold off
title('Error')
xlabel('x');ylabel('Error')
legend('RK4 Error','ABM4 Error','Zero reference','Location','Best')
fprintf('RK1\n')
disp(y1RK(:,1:4))
fprintf('ABM1\n')
disp(y1ABM(:,1:4))

%% Setup 2: Spring mass damper system
% Tuning parameters
c = 1; % Damping constant
k = 1; % Spring constant
m = 1; % Mass 
h2 = 0.01;
%h2 = sqrt(abs(-k/m-(c/m)^2/4));
x0_2 = 0;
xf_2 = 10;
y0_2 = [2;0];
% EoM
eqn2 = @(t,x) [x(2);-c/m*x(2)-k/m*x(1)];
% Numerical Solution
tspan_2 = [x0_2, xf_2];
sol_2RK = RK4(eqn2,tspan_2,y0_2,h2);
sol_2ABM = ABM4(eqn2,tspan_2,y0_2,h2);
% Extract solutions
x2RK = sol_2RK.x;
y2RK = sol_2RK.y(1,:);
x2ABM = sol_2ABM.x;
y2ABM = sol_2ABM.y(1,:);
% Exact solution
omega = sqrt(abs(-k/m-(c/m)^2/4));
y2_exactRK = exp(-c/(2*m)*x2RK).*(y0_2(1)*cos(omega*x2RK)+...
    2*y0_2(2)*m/(c*omega)*sin(omega*x2RK));
y2_exactABM = exp(-c/(2*m)*x2ABM).*(y0_2(1)*cos(omega*x2ABM)+...
    2*y0_2(2)*m/(c*omega)*sin(omega*x2ABM));
% Error
er2RK = -y2_exactRK + y2RK;
er2ABM = -y2_exactABM + y2ABM;
% Plotting
figure
subplot(2,1,1)
plot(x2RK,y2RK,x2ABM,y2ABM,'-.',x2RK,y2_exactRK,'--','linewidth',1.5)
title('RK4 vs ABM4 vs Exact solution')
xlabel('x');ylabel('y')
legend('RK4','ABM4','Exact solution','Location','best')
subplot(2,1,2)
plot(x2RK,er2RK,x2ABM,er2ABM,'--','linewidth',2.5)
hold on
plot(tspan_2,[0 0],'--k')
hold off
title('Error')
xlabel('x');ylabel('Error')
legend('RK4 Error','ABM4 Error','Zero reference','Location','Best')
fprintf('RK2\n')
disp(sol_2RK.y(:,1:4))
fprintf('ABM2\n')
disp(sol_2ABM.y(:,1:4))

%% Setup 3 Double pendulum
% Tuning parameters
x0_3 = 0;
xf_3 = 7;
y0_3 = [60*pi/180;60*pi/180;2;8];
P.g = 9.81; % Gravity
P.m = [1 1]; % Mass vector
P.l = [1;1]; % Massless rod length vector
h3 = .01;
% Numerical Solution
tspan_3 = [x0_3, xf_3];
sol_3RK = RK4(@(t,x)eqn(t,x,P),tspan_3,y0_3,h3);
sol_3ABM = ABM4(@(t,x)eqn(t,x,P),tspan_3,y0_3,h3);
sol_345 = ode45(@(t,x)eqn(t,x,P),tspan_3,y0_3);
% Extract solutions
x3RK = sol_3RK.x;
y3RK1 = sol_3RK.y(1,:);
y3RK2 = sol_3RK.y(2,:);
y3RK1d = sol_3RK.y(3,:);
y3RK2d = sol_3RK.y(4,:);
x3ABM = sol_3ABM.x;
y3ABM1 = sol_3ABM.y(1,:);
y3ABM2 = sol_3ABM.y(2,:);
y3ABM1d = sol_3ABM.y(3,:);
y3ABM2d = sol_3ABM.y(4,:);
x345 = sol_345.x;
y3451 = sol_345.y(1,:);
y3452 = sol_345.y(2,:);
y3451d = sol_345.y(3,:);
y3452d = sol_345.y(4,:);
% Plotting
figure
plot(y3RK1,y3RK1d,y3ABM1,y3ABM1d,'-.',y3451,y3451d,'--','linewidth',1.5)
title('RK4 vs ABM4 vs ODE45 for \theta_1 vs \theta''_1')
xlabel('\theta');ylabel('\theta''')
legend('RK4','ABM4','ODE45','Location','best')
figure
plot(y3RK2,y3RK2d,y3ABM2,y3ABM2d,'-.',y3452,y3452d,'--','linewidth',1.5)
title('RK4 vs ABM4 vs ODE45 for \theta_2 vs \theta''_2')
xlabel('\theta');ylabel('\theta''')
legend('RK4','ABM4','ODE45','Location','best')