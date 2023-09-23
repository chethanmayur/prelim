clear all;
close all;
clc;


PO1 = 101353;                   % inlet stagnation pressure
TO1 = 288.15;                   % inlet stagnation temperature
m = 0.8;                        % mass flow rate
N = 90790;                      % RPM
Cp = 1005;                      % Specific heat
k = 1.4;                        % gamma
R = 287;                        % gas constant
D2 = 0.1;                       % exit diameter
ZB = 18;                         % no of blades
alpha1 = 0;                     % absolute inlet velocity angle
Rh1 = 0.015;                     % hub radius 
B1 = 0.04;                      % boundary layer blackage
lambda2 = 2;                    % exit swirl parameter
eta = 0.78;                      % stage efficiency
beta2b = -40;                   % exit relative dirxn degree
AR = 2.5;                       % area ratio
pr = 4;                         % pressure ratio


Cm1i = [100:0.5:300]';          % inlet absolute meridonial velocity range
TOi = [293*ones(length(Cm1i))]; % inlet stagnation temperature range


Ctheta1i = Cm1i*tand(alpha1);   % Inlet absolute tangential velocity [degrees]

C1i = (Ctheta1i.^2 + Cm1i.^2).^0.5; % inlet absolute velocity

T1i = TO1 - (C1i.^2/(2*Cp));    % Inlet temperature

P1i = PO1 * (T1i/TO1).^(k/(k-1)); % Inlet pressure

sound1i = (k*R*T1i).^0.5;
M1i = T1i./sound1i;                  % Inlet mach no

rho1i = P1i./(R*T1i);            % Inlet density

A1i = m./(rho1i.*Cm1i*(1-B1));       % inlet area

R1i = ((A1i/(pi())) + Rh1^2).^0.5;    % inlet tip radius

U1i = 2*pi()*R1i*N/60;              % blade tip speed

W1ti = (Cm1i.^2 +((U1i.^2)-(Ctheta1i.^2))).^0.5; % inlet relative velocity

x = find(W1ti==min(W1ti));           % find pos of the minimum relative vel


W1 = W1ti(x);
C1 = C1i(x)
T1 = T1i(x);
M1 = M1i(x);
P1 = P1i(x);
rho1 = rho1i(x);
A1 = A1i(x);
R1 = R1i(x);
U1 = U1i(x);
Cm1 = Cm1i(x);
Ctheta1 = Ctheta1i(x);

beta1 = atand((U1-Ctheta1)/Cm1);         % inlet relative vel angle

%sigma = 1 - ((cosd(beta2b))^0.5/(ZB)^0.7); % slip factor
sigma = 1-(pi()*0.63/ZB)                    % slip factor

mu = sigma*lambda2/(lambda2-tand(beta2b));  % work input coefficient

hx = ((k*R*TO1)/(k-1))*(pr^((k-1)/k) -1);   % specific enthalpy

wx = hx/eta;            % specific work

TO2 = TO1 + (wx*(k-1)/(k*R));            % stagnation exit temp

PO2 = PO1 * (((k-1)*wx*eta/(k*R*TO1))+1)^(k/(k-1))      % stagnation exit press

U2 = D2*pi()*N/60;                   % exit blade speed

% U2 = ((U1*Ctheta1+wx)/mu)^0.5 #Exit Blade Speed 

Ctheta2 = mu*U2;                            % Absolute tangential exit velocity

Cm2 = Ctheta2/lambda2;                        % Absolute meridional exit velocity

T2 = TO2 - (k-1)*(Ctheta2^2+Cm2^2)/(2*k*R)       %Exit temperature

P2 = PO2 / ((TO2/T2)^(k/(k-1)));                %exit pressure

rho2 = P2/(T2*R);             % Exit density [kg/m^3]

A2 = m/(rho2*Cm2);            % Exit area

b2 = A2/(pi()*D2)          % Depth of impeller exit

C2 = (Ctheta2^2+Cm2^2)^0.5          % Absolute exit velocity


etad = 0.85; %Estimated diffuser efficiency
CpDi = 1-(AR^-2); %Ideal pressure recovery coefficient
CpD = etad*CpDi; %Pressure recovery coefficient
P3 = P2 + CpD*(PO2-P2); %Diffuser exit static pressure [Pa]
C3 = C2/AR %Diffuser exit absolute velocity [m/s]
PO3 = P3 +0.5*rho2*C3^2 %Diffuser exit stagnation pressure [Pa]


etaiterate=((PO3/PO1)^((k-1)/k)-1)/((TO2/TO1)-1) %Iterative stage efficiency [-]
Prest=((etaiterate*U2^2*mu)/(Cp*T1)+1)^(k/(k-1)) %Estimate of the pressure ratio

work = wx*m                                        % work done