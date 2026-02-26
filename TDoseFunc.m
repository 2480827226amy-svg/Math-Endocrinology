%% Dosing function iteration 3 
% Incorporate a weekly subQ dose into the EG model
% Goal: add in a thirteenth function, Adipose T, through which the T
% dose enters in and then slowly diffuses out into Serum T.

% Unfinished as of 5/20/25

%Initial Conditions
RFSH0 = 116.82;
SFSH0 = 142.5;
RLH0 = 250.35;
SLH0= 25.34;
dphi0 = 0.50185;
domega0 = 9.7509;
dlambda0 = 4.102;
S0= 0.050498;
SerumT0 = 273.67;
InterT0 = 0.003999;
E20=56.387;
P40 = 0.468;
AdpT0 = 0;

y0 = [RFSH0, SFSH0, RLH0, SLH0, dphi0, domega0, dlambda0, S0, SerumT0, InterT0, AdpT0, E20, P40];

for  ts = 1:50
Tdose = 8*10^7; %in ng 

tspan1 = [ts ts+1];
opts = odeset('AbsTol', 1e-6,'RelTol', 1e-8,'MaxStep', 0.1);
[t, y] = ode15s(@(t,y) f(t,y), tspan1, y0);

%FOR A WEEKLY DOSE 
if mod(ts, 7) == 0
    y0 = [(y(end,1:10)) (y(end,11) + Tdose) (y(end,12:13))];
    AdpT0 = y(end,11)+ Tdose;
else 
   y0 = [(y(end,1:10)) (y(end,11)) (y(end,12:13))];
   AdpT0 = y(end,11);
end

figure(1)
plot(t,y(:,9),"LineWidth",1,"Color","black")

hold on
end
hold off

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
AdpT = y(11);
E2 = y(12);
P4 = y(13);

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
SerumT0 = 273.67;
Tdose = 8*10^7;
AdpT0 = y(11) + Tdose; 
%Rate of diffusion from fat to bloodstream
r1 = 1;
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
dydx(9) = t0 + 2 - deltaT * SerumT + (((t1 * G1) * (F1 + cTF2 * F2)) + t2 * G1 * G2 * F1) * (dphi + tau1 * domega + tau2 * S * dlambda + tau3 * (1 - ((dphi + domega + dlambda) / psi)));
% ASSUMING G1 and G2 to be basal levels of insulin so G1=G2=1
%Inter T
dydx(10) = (tg1 * G1 * G2 * F1) - (((tg2 * SFSH) / (h3 + SFSH)) * (InterT));
%AdpT
dydx(11) = -2;
% Serum E2
dydx(12) = e0 - (deltaE * E2) + ((tg2 * SFSH) / (h3 + SFSH)) * (InterT) * (dphi + (eta * dlambda * S));
% Serum P4
dydx(13) = (-deltaP) * P4 + ((p * SLH) / (SLH + hp)) * (S * dlambda);
dydx = dydx';
end