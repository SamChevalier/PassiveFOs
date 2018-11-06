%% 3 Element System Simulation
%  1 Bus: Infinite
%
%
%    Inf_______________Inf
%         |    |    |
%        Gen   Z    PQ
%         |    |    |
%         V    V    V
%
clear variables;

%% Simulation Parameters
j = sqrt(-1);

% Simulation Length
dt       = 0.01;     % Time Step
tf       = 60;       % Simulation Length
t_vec    = 0:dt:tf;
w        = 2*pi*1;   % Oscillation Frequnecy

%% ------------ Test 1 ------------ %%

% System Parameters
Vm_inf = 1; % Infinite bus voltage magnitude
Th_inf = 0; % Infinite bus voltage angle

% Generator Data (Connection 1)
M1  = 2.5;
D1  = 5;
E1  = 1;
Xd1 = 0.1;

% Impedance Data (Connection 2)
R2 = 10;
X2 = 1;
G2 = real(1/(R2 + j*X2));

% PQ Load Data (Connection 3)
P3 = 0.5;
Q3 = 0.25;

% (Vr perturbations) lead (Vi perturbations)
M_Vr = 0.1;   % Magnitude of Perturbation
M_Vi = 0.1;

P_Vr = pi/10;  % Phase of Perturbation
P_Vi = -pi/10;

% Impedance and PQ Loads are algebraic: compute explicitly
V_inf = Vm_inf*exp(j*Th_inf) + M_Vr*sin(w*t_vec + P_Vr) + j*M_Vi*sin(w*t_vec + P_Vi);
dVr   = M_Vr*sin(w*t_vec + P_Vr);
dVi   = M_Vi*sin(w*t_vec + P_Vi);

% Derivative (analytical)
d_dVr   = w*M_Vr*cos(w*t_vec + P_Vr);
d_dVi   = w*M_Vi*cos(w*t_vec + P_Vi);

% Impedance Load
Iz     = V_inf/(R2+j*X2);
dIz_r2 = real(Iz) - real(Vm_inf*exp(j*Th_inf)/(R2+j*X2));
dIz_i2 = imag(Iz) - imag(Vm_inf*exp(j*Th_inf)/(R2+j*X2));

% PQ Load
Ipq     = conj((P3+j*Q3)./(V_inf));
dIpq_r3 = real(Ipq) - real((P3+j*Q3)./(Vm_inf*exp(j*Th_inf)));
dIpq_i3 = imag(Ipq) - imag((P3+j*Q3)./(Vm_inf*exp(j*Th_inf)));

% Generator Response
syms w1(t) d1(t) Ti(t) Vi(t)

% ODE Variables
ODEvars = [w1(t) d1(t)];

% Electrical Powers
Pe1 = (E1*Vi(t)/Xd1)*sin(d1(t) - Ti(t));

% Del and Pm IC
IC.d1 = Th_inf + pi/4;
Pm1   = (E1*Vm_inf/Xd1)*sin(IC.d1 - Th_inf);

% ODEs
ODEs = [diff(w1(t))    == (Pm1 - Pe1 - D1*w1(t))/M1, ...
        diff(d1(t))    == w1(t)];

% Set Up System
ff  = daeFunction(ODEs, ODEvars, Ti(t), Vi(t));

% Voltage Magnitude and Phase Ocillation
Ti = @(t) angle(Vm_inf*exp(j*Th_inf) + M_Vr*sin(w*t + P_Vr) + j*M_Vi*sin(w*t + P_Vi));
Vi = @(t) abs(  Vm_inf*exp(j*Th_inf) + M_Vr*sin(w*t + P_Vr) + j*M_Vi*sin(w*t + P_Vi));

% Redefine System of Equations
FF   = @(t, Y, YP) ff(t, Y, YP, Ti(t), Vi(t));

% Assign the initial conditions
yp0est = zeros(2,1);
y0est  = [0;
          IC.d1];

% Determine True Initial Conditions
[y0, yp0] = decic(FF, 0, y0est, [], yp0est, []);

% Simulate the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan  = 0:dt:tf;
[t_out,y_out]  = ode15i(FF,tspan,y0,yp0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r.w1  = y_out(:,1);
r.d1  = y_out(:,2);

% Get Currents (Flow into Generator)
Ig   = (V_inf.'-E1*exp(j*r.d1))/(j*Xd1);
Ig_0 = (Vm_inf*exp(j*Th_inf)-E1*exp(j*IC.d1))/(j*Xd1);
dIg_r1 = real(Ig) - real(Ig_0);
dIg_i1 = imag(Ig) - imag(Ig_0);

% Compute DEF Integrals
W1_gen = cumtrapz((dIg_r1.'.*d_dVi - dIg_i1.'.*d_dVr)*dt);
W2_z   = cumtrapz((dIz_r2.*d_dVi - dIz_i2.*d_dVr)*dt);
W3_pq  = cumtrapz((dIpq_r3.*d_dVi - d_dVr.*dIpq_i3)*dt);

% Generate Plots
clf
c1  = [0         0.4470    0.7410];
c1s = [0         0.4470    0.7410   0.8];
c2  = [0.8500    0.3250    0.0980];
c2s = [0.8500    0.3250    0.0980   0.2];
b   = [0 0 0];
bs  = [0 0 0 0.8];

% Data Drive Plots
hold on
plot(t_vec,W1_gen,'Linewidth',0.8,'color',c1s);
hold on
plot(t_vec,t_vec*sin(P_Vr - P_Vi)*M_Vr*M_Vi*(w^2)*G2*2/(w*2),'Linewidth',6,'color',c2s);
plot(t_vec,W2_z,'Linewidth',1,'color',c2);
plot(t_vec,W3_pq,'Linewidth',0.8,'color',bs);

xlabel('$\rm{Time\:(sec)}$','Interpreter','latex','FontSize',15)
ylabel('${\rm Dissipating \;\, Energy}$','Interpreter','latex','FontSize',15)
legend({'$E^{\star}_{g}$';'$E^{\star}_{z,\;{\rm prediction}}$';'$E^{\star}_{z}$';'$E^{\star}_{p}$'},'Location','northwest','Interpreter','latex','FontSize',15,'Box','off');
set(gca,'FontName','Times','FontSize',15)
set(gcf,'Units','inches','Position',[0 0 9 3])
xlim([0 60])
tightfig(gcf); % Now "File" => "Export Setup" => "Expand axes to fill figure"


%% ------------ Test 2 ------------ %%

% System Parameters
Vm_inf = 1; % Infinite bus voltage magnitude
Th_inf = 0; % Infinite bus voltage angle

% Generator Data (Connection 1)
M1  = 2.5;
D1  = -0.3;
E1  = 1;
Xd1 = 0.1;

% Impedance Data (Connection 2)
R2 = 10;
X2 = 1;
G2 = real(1/(R2 + j*X2));

% PQ Load Data (Connection 3)
P3 = 0.5;
Q3 = 0.25;

% (Vr perturbations) lead (Vi perturbations)
M_Vr = 0.1;   % Magnitude of Perturbation
M_Vi = 0.1;
P_Vr = -pi/10;  % Phase of Perturbation
P_Vi = pi/10;

% Impedance and PQ Loads are algebraic: compute explicitly
V_inf = Vm_inf*exp(j*Th_inf) + M_Vr*sin(w*t_vec + P_Vr) + j*M_Vi*sin(w*t_vec + P_Vi);
dVr   = M_Vr*sin(w*t_vec + P_Vr);
dVi   = M_Vi*sin(w*t_vec + P_Vi);

% Derivative (analytical)
d_dVr   = w*M_Vr*cos(w*t_vec + P_Vr);
d_dVi   = w*M_Vi*cos(w*t_vec + P_Vi);

% Impedance Load
Iz     = V_inf/(R2+j*X2);
dIz_r2 = real(Iz) - real(Vm_inf*exp(j*Th_inf)/(R2+j*X2));
dIz_i2 = imag(Iz) - imag(Vm_inf*exp(j*Th_inf)/(R2+j*X2));

% PQ Load
Ipq     = conj((P3+j*Q3)./(V_inf));
dIpq_r3 = real(Ipq) - real((P3+j*Q3)./(Vm_inf*exp(j*Th_inf)));
dIpq_i3 = imag(Ipq) - imag((P3+j*Q3)./(Vm_inf*exp(j*Th_inf)));

% Generator Response
syms w1(t) d1(t) Ti(t) Vi(t)

% ODE Variables
ODEvars = [w1(t) d1(t)];

% Electrical Powers
Pe1 = (E1*Vi(t)/Xd1)*sin(d1(t) - Ti(t));

% Del and Pm IC
IC.d1 = Th_inf + pi/4;
Pm1   = (E1*Vm_inf/Xd1)*sin(IC.d1 - Th_inf);

% ODEs
ODEs = [diff(w1(t))    == (Pm1 - Pe1 - D1*w1(t))/M1, ...
        diff(d1(t))    == w1(t)];

% Set Up System
ff  = daeFunction(ODEs, ODEvars, Ti(t), Vi(t));

% Voltage Magnitude and Phase Ocillation
Ti = @(t) angle(Vm_inf*exp(j*Th_inf) + M_Vr*sin(w*t + P_Vr) + j*M_Vi*sin(w*t + P_Vi));
Vi = @(t) abs(  Vm_inf*exp(j*Th_inf) + M_Vr*sin(w*t + P_Vr) + j*M_Vi*sin(w*t + P_Vi));

% Redefine System of Equations
FF   = @(t, Y, YP) ff(t, Y, YP, Ti(t), Vi(t));

% Assign the initial conditions
yp0est = zeros(2,1);
y0est  = [0;
          IC.d1];

% Determine True Initial Conditions
[y0, yp0] = decic(FF, 0, y0est, [], yp0est, []);

% Simulate the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan  = 0:dt:tf;
[t_out,y_out]  = ode15i(FF,tspan,y0,yp0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r.w1  = y_out(:,1);
r.d1  = y_out(:,2);

% Get Currents (Flow into Generator)
Ig   = (V_inf.'-E1*exp(j*r.d1))/(j*Xd1);
Ig_0 = (Vm_inf*exp(j*Th_inf)-E1*exp(j*IC.d1))/(j*Xd1);
dIg_r1 = real(Ig) - real(Ig_0);
dIg_i1 = imag(Ig) - imag(Ig_0);

% Compute DEF Integrals
W1_gen = cumtrapz((dIg_r1.'.*d_dVi - dIg_i1.'.*d_dVr)*dt);
W2_z   = cumtrapz((dIz_r2.*d_dVi - d_dVr.*dIz_i2)*dt);
W3_pq  = cumtrapz((dIpq_r3.*d_dVi - d_dVr.*dIpq_i3)*dt);

% Generate Plots
clf
c1  = [0         0.4470    0.7410];
c1s = [0         0.4470    0.7410   0.8];
c2  = [0.8500    0.3250    0.0980];
c2s = [0.8500    0.3250    0.0980   0.2];
b   = [0 0 0];
bs  = [0 0 0 0.8];

% Data Drive Plots
hold on
plot(t_vec,W1_gen,'Linewidth',0.8,'color',c1s);
hold on
plot(t_vec,t_vec*sin(P_Vr - P_Vi)*M_Vr*M_Vi*(w^2)*G2*2/(w*2),'Linewidth',6,'color',c2s);
plot(t_vec,W2_z,'Linewidth',1,'color',c2);
plot(t_vec,W3_pq,'Linewidth',0.8,'color',bs);

xlabel('$\rm{Time\:(sec)}$','Interpreter','latex','FontSize',15)
ylabel('${\rm Dissipating \;\, Energy}$','Interpreter','latex','FontSize',15)
legend({'$E^{\star}_{g}$';'$E^{\star}_{z,\;{\rm prediction}}$';'$E^{\star}_{z}$';'$E^{\star}_{p}$'},'Location','southwest','Interpreter','latex','FontSize',15,'Box','off');
set(gca,'FontName','Times','FontSize',15)
set(gcf,'Units','inches','Position',[0 0 9 3])
ylim([-0.3 0.1])
tightfig(gcf); % Now "File" => "Export Setup" => "Expand axes to fill figure"