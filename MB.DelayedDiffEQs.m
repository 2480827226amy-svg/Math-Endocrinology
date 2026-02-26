%% Exploration of a delayed model
% Drs Mukhopadhyay and Bhattacharyya
% Using delayed differential equations will account for the time 
% it takes for the T to move through the body.

% Randomly assigned time values
tau1 = 3;
tau2 = 1;
tau0 = 1;
lags = [tau1 (tau0+tau2)];  % Lags corresponding to L and T
tspan = [0 5];

sol = dde23(@model2, lags, @history, tspan);

plot(sol.x,sol.y,'-o')
title('A delayed mathematical model for testosterone secretion with feedback control mechanism')
subtitle('by Drs Mukhopadhyay and Bhattacharyya')
xlabel('t days')
ylabel('Solutions')
legend('LHRT','LH','T')

function dydt = model2(t,y,Z);
d1 = 0.4;
d2 = 0.4;
d3 = 0.4;
r2 = .8;
r3 = .8;
a = 0.0005;
Tdose = 10;
%randomly assigned parameters 

ylag1 = Z(:,1);
ylag2 = Z(:,2);

dydt = [%R - LHRT secretion:
-a*y(3) - d1*y(1);
%L - LH secretion:
r2* ylag1(1) - d2*y(2);
%T - T secretion:
r3*ylag2(2) - d3*y(3) + Tdose];

end

function s = history(t);
% Define a history function that returns initial conditions
s = [20,10,3];
end
