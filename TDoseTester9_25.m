%% Our final model for incorporating exogenous T into Dr. Graham's model
% Note that the T is incorporated by adding it into the initial conditions

%Initial Conditions
%All of these are the same as those in EG paper
RFSH0 = 116.82;
SFSH0 = 142.5;
RLH0 = 250.35;
SLH0= 25.34;
dphi0 = 0.50185;
domega0 = 9.7509;
dlambda0 = 4.102;
S0= 0.050498;
SerumT0= 273.67;
InterT0 = 0.003999;
E20 = 56.387;
P40 = 0.468;

%Sets up a vector of all of the initial conditions
y0 = [RFSH0, SFSH0, RLH0, SLH0, dphi0, domega0, dlambda0, S0, SerumT0, InterT0, E20, P40];

%Solves the system of ODEs 
% Note the SerumT0 at the end, updates the value with a new dose each day

time = 365;

figure
%PLOT 1: NO DOSE
for  ts = 0:time
hold on
Tdose = 273; 
tspan1 = [ts ts+1];
opts = odeset('AbsTol', 1e-6,'RelTol', 1e-8,'MaxStep', 0.1);
[t, y] = ode15s(@(t,y) f(t,y), tspan1, y0);
y0 = [(y(end,1:8)) (y(end,9)+Tdose) (y(end,10:12))];
SerumT0 = y(end,9)+Tdose;

plot1 = plot(t,y(:,4),"DisplayName","No Dose","Color","#6495ED","LineStyle","-","LineWidth",1.5);

% %Plotting the solutions: focused on LH and P4 specifically
% % LH graph
% subplot(2,1,1)
% plot(t,y(:,4),"Color","#6495ED","LineStyle","-","LineWidth",1.5)
% xlabel('Time (days)')
% ylabel('LH (\mu g/L)')
% fontsize(.6,"centimeters")
% hold on
% % P4 graph
% subplot(2,1,2)
% plot(t,y(:,12),"Color","#6495ED","LineStyle","-","LineWidth",1.5)
% xlabel ('Time (days)')
% ylabel('P4 (ng/mL)')
% fontsize(.6,"centimeters")

end


%PLOT 2: LOLO DOSE
for  ts = 0:time
hold on
Tdose = 273*10^2; %MIDDLE
tspan1 = [ts ts+1];
opts = odeset('AbsTol', 1e-6,'RelTol', 1e-8,'MaxStep', 0.1);
[t, y] = ode15s(@(t,y) f(t,y), tspan1, y0);
y0 = [(y(end,1:8)) (y(end,9)+Tdose) (y(end,10:12))];
SerumT0 = y(end,9)+Tdose;

%Plotting LH
plot2 = plot(t,y(:,4),"DisplayName","LoLo Dose","Color","r","LineStyle","-","LineWidth",1.5);
end


%PLOT 3: LOW DOSE
for  ts = 0:time
hold on
Tdose = 1.79*10^6; %LOW DOSE IN PAPER 
tspan1 = [ts ts+1];
opts = odeset('AbsTol', 1e-6,'RelTol', 1e-8,'MaxStep', 0.1);
[t, y] = ode15s(@(t,y) f(t,y), tspan1, y0);
y0 = [(y(end,1:8)) (y(end,9)+Tdose) (y(end,10:12))];
SerumT0 = y(end,9)+Tdose;

%Plotting LH
plot3 = plot(t,y(:,4),"DisplayName","Low Dose","Color","g","LineStyle","-","LineWidth",1.5);
end


%PLOT 4: HIGH DOSE
for  ts = 0:time
hold on
Tdose = 3.57*10^6; %HIGH DOSE IN PAPER 
tspan1 = [ts ts+1];
opts = odeset('AbsTol', 1e-6,'RelTol', 1e-8,'MaxStep', 0.1);
[t, y] = ode15s(@(t,y) f(t,y), tspan1, y0);
y0 = [(y(end,1:8)) (y(end,9)+Tdose) (y(end,10:12))];
SerumT0 = y(end,9)+Tdose;

%Plotting LH
plot4 = plot(t,y(:,4),"DisplayName","High Dose","Color","k","LineStyle","-","LineWidth",1.5);
end

legend([plot1 plot2 plot3 plot4])


% Setting up the system of ODEs
function  dydx = f(~, y)
RFSH = y(1);
SFSH = y(2);
RLH = y(3);
SLH = y(4);
dphi = y(5);
domega = y(6);
dlambda = y(7);
S = y(8);
SerumT = y(9);
InterT = y(10);
E2 = y(11);
P4 = y(12);
% All the parameters, grouped
% Group A
cFE = 0.0022729;
cFI = 1.9488;
cFP = 60.428;
cLE = 0.0010404;
cLP = 0.0099415;
cLT = 0.0095942;
deltaF = 8.21;
deltaL = 14;
kF = 2.5412;
KFI = 107.01;
KiLP = 0.34952;
KLT = 420;
KmL = 183.56;
kL = 0.74567;
n = 8;
V = 2.5;
v0L = 1051.7;
v1L = 34838;
vF = 3236.6;
% Group B
cphiF = 0.01127;
cphiT = 0.19878;
deltaS = 0.74702;
f0 = 0.0025112;
f1 = 4.3764;
f2 = 27.812;
h1 = 590.32;
h2 = 1815.3;
hp = 20.764;
hs = 12.329;
l = 0.49017;
m = 4;
shat = 2.378;
w = 0.23173;
% Group C
cTF2 = 123.8136;
deltaE = 1.1;
deltaP = 0.5;
deltaT = 5.5;
e0 = 44.512;
eta = 1.1087;
h3 = 17.796;
k1 = 1.09;
k2 = 22.28645;
k3 = 113.9188;
p = 0.3734;
t0 = 741.68;
t1 = 0.57088;
t2 = 1.3481;
tau1 = 5.3989;
tau2 = 0;
tau3 = 430.91;
tg1 = 6.6548;
tg2 = 186.27;
psi = 2004.3;
G1 = 1;
G2 = 1;
D = k1 * (SLH^2) + k2 * SLH + k3;
F1 = (SLH^2) / D;
F2 = SLH / D;
Tdose = 2.9*10^6;
SerumT0= y(9)+Tdose; 
%273.67; (I don't know what this number is)

%% The differential equations themselves:
%% Subsystem 1
% Releasable FSH, FSHp
dydx(1) = (vF / (1 + (cFI * ((S * dlambda) / (KFI + (S * dlambda)))))) - kF * ((1 + (cFP * P4)) / (1 + cFE * (E2^2))) * RFSH;
% Serum FSH
dydx(2) = (1 / V) * kF * ((1 + cFP * P4) / (1 + cFE * (E2)^2)) * RFSH - (deltaF * SFSH);
% Releasable LH, LHp
dydx(3) = (v0L * (SerumT / (KLT + SerumT)) + v1L * ((E2)^n / ((KmL)^n + (E2)^n))) * (1 / (1 + (P4 / (KiLP * (1 + cLT * SerumT))))) - kL * ((1 + cLP * P4) / (1 + cLE * E2)) * RLH;
% Serum LH
dydx(4) = (1 / V) * kL * ((1 + cLP * P4) / (1 + cLE * E2)) * RLH - deltaL * SLH;
%% Subsystem 2
% Follicular Phase, phi
dydx(5) = f0 * (SerumT / (SerumT0)) + dphi * (((f1 * SFSH^2) / ((h1 / (1 + cphiT * (SerumT / (SerumT0))))^2 + SFSH^2)) - ((f2 * SLH^2) / ((h2 / (1 + cphiF * SFSH))^2 + SLH^2)));
% Ovulatory phase, sigma
dydx(6) = ((f2 * SLH^2) / ((h2 / (1 + cphiF * SFSH))^2 + SLH^2)) * dphi - w * S * domega;
% Luteal Phase, lambda
dydx(7) = (w * S * domega) - l * (1 - S) * dlambda;
% LH Support
dydx(8) = ((shat * SLH^m) / (hs^m + SLH^m)) * (1 - S) - deltaS * S;
%% Subsystem 5 combined
% Serum T, T
dydx(9) = t0 - deltaT * SerumT + (((t1 * G1) * (F1 + cTF2 * F2)) + t2 * G1 * G2 * F1) * (dphi + tau1 * domega + tau2 * S * dlambda + tau3 * (1 - ((dphi + domega + dlambda) / psi)));
% Intermediate T, T-gamma (replacing granulosa t from original subsystem 3)
% ASSUMING G1 and G2 to be basal levels of insulin so G1=G2=1
%Inter T
dydx(10) = (tg1 * G1 * G2 * F1) - (((tg2 * SFSH) / (h3 + SFSH)) * (InterT));
% Serum E2
dydx(11) = e0 - (deltaE * E2) + ((tg2 * SFSH) / (h3 + SFSH)) * (InterT) * (dphi + (eta * dlambda * S));
% Serum P4
dydx(12) = (-deltaP) * P4 + ((p * SLH) / (SLH + hp)) * (S * dlambda);
dydx = dydx';
end