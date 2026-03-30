%% 

close all
clc
clear

%% constants
global g G m_earth r_earth mu h_orbit v_earth vehicle t_seperation

g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m
mu=G*m_earth;
omega_earth=7.2921159e-5; %units=rads/s

%% Loop

vehicle=Atlas_V_501();
altitude=vehicle.performance.altitude;
interval=vehicle.performance.interval;
payload=vehicle.performance.payload;
h_vec = altitude(1):interval:altitude(end);     % meters
lat_launch=deg2rad(28); %cape canaveral
stage=vehicle.stage;

N = numel(h_vec);
m_pld_max_vec = nan(N,1);       % payload result (kg)
h_apo_vec     = nan(N,1);       % achieved apogee (m)
h_per_vec     = nan(N,1);       % achieved perigee (m)
e_vec         = nan(N,1);       % eccentricity
data_all(1:N) = struct('t',[],'h',[],'v',[],'m',[],'x',[],'theta',[]);

for k = 1:N
    
%% Setup

%specify mission 
h_orbit = h_vec(k);

azimuth=deg2rad(90); %due east 

v_earth=omega_earth*r_earth*cos(lat_launch); %pure eastward velocity
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

t_seperation=15;

burn_time_max=stage(2).m_prop/(stage(2).thrust/(stage(2).isp*g));
burn_time_guess=300+((burn_time_max-300)/2);
pld_guess=stage(1).pld_guess+((stage(2).pld_guess-stage(1).pld_guess)/2);

%% Optimisation

%scale variables to around 1 
sol_scale = [1000, 0.01, 10, 10, 100, 0.1];

%initial guesses
sol_0=   [ 
    pld_guess; %m_pld
    deg2rad(2); %theta_kick
    15; %t_pitchover
    15; %t_pitch
    burn_time_guess; %burn_time_2
    deg2rad(2) %theta_SECO
]./sol_scale';

lowbound=[ 
    stage(1).pld_guess; %m_pld
    deg2rad(1); %theta_kick
    5; %t_pitchover
    10; %t_pitch
    300; %burn_time_2
    deg2rad(0) %theta_SECO
]./sol_scale';

upbound= [ 
    stage(2).pld_guess; %m_pld
    deg2rad(5); %theta_kick
    20;  %t_pitchover
    20; %t_pitch
    burn_time_max;  %burn_time_2
    deg2rad(10)  %theta_SECO
]./sol_scale'; 

objective=@(sol) -sol(1)*sol_scale(1);

nonlcon=@(sol) constraints(sol,sol_scale,h_orbit);

options=optimoptions('fmincon', ...
    'Algorithm','interior-point', ...
    'FiniteDifferenceType','central' , ...
    'Display','iter', ...
    'FiniteDifferenceStepSize',1e-3, ...
    'ConstraintTolerance',1e-4, ...
    'OptimalityTolerance', 1e-4, ...
    'MaxIterations', 2000 ...
    );

sol__optimal_scaled=fmincon(objective,sol_0,[],[],[],[],lowbound,upbound,nonlcon,options);

sol_optimal=sol__optimal_scaled.*sol_scale';



%% Results 

data=sim(sol_optimal);
[h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit,sol_optimal);

m_pld_max = sol_optimal(1) + (m_prop_remaining - m_prop_2b);

% store
m_pld_max_vec(k) = m_pld_max;
h_apo_vec(k)     = h_apo;
h_per_vec(k)     = h_per;
e_vec(k)         = e;
data_all(k)      =data;


end

accuracy=100.*((m_pld_max_vec-payload)./payload);

figure;
plot(h_vec'./1000, (h_apo_vec - h_vec')/1000, '-+'); 
hold on
plot(h_vec'./1000, (h_per_vec - 120e3)/1000, '-o');
grid
xlabel('Target orbit altitude (km)')
ylabel('Constraint residual (km)')
legend('Apogee error', 'Perigee margin')
hold off

figure;
plot(h_vec/1000, m_pld_max_vec, '-o');
grid
hold on
plot(altitude/1000,payload);
legend('Simulated','Actual')
xlabel('Target orbit altitude (km)')
ylabel('Max payload (kg)')
hold off


figure; hold on; grid on

for k = 1:N
    plot(data_all(k).t, data_all(k).theta);
end

xlabel('Time (s)')
ylabel('Theta (rad)')
hold off

figure;
bar(h_vec'./1000,accuracy)
grid
xlabel('Target orbit altitude (km)')
ylabel('Accuracy (%)')
hold on
yline(10,'r');
yline(-10,'r');
hold off

%% Script functions

function [c,ceq] = constraints(sol,sol_scale,h_orbit)

sol_unscaled=sol.*sol_scale';

data =sim(sol_unscaled);

[h_apo, h_per, ~, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit, sol_unscaled);

h_min=120e3;
h_cut_min = 150e3;         
h_SECO=data.h(end);

    % equality constraints: hit target altitude, mass
    ceq =  [(h_apo - h_orbit)/1000;
           ];

    c = [(m_prop_2b-m_prop_remaining)/10;
         (h_min-h_per)/1000;
         (h_cut_min - h_SECO)/1000
    ];
end 

function[rho_atm,P_atm] = isa_model(h)

R=287; %units=J/kg*K
g=9.80665; %units=m/s^2

h = max(h,0);

%ISA model
if h<11e3 %troposphere
    T_atm =288.15-(h *6.5e-3);
    P_atm =101325*(T_atm /288.15)^(g/(R*(6.5e-3)));
elseif h >=11e3 && h <=20e3 %tropopause    
    T_atm =216.65;
    P_atm =22632*exp((-1*g*(h -11000)/(R*T_atm )));
elseif h >20e3 && h <32e3 %stratosphere (lower)
    T_atm =216.65+((h -20e3)*1e-3);
    P_atm =5474*(T_atm /216.65)^(g/(R*-1e-3));
elseif h >=32e3 && h <47e3 %stratosphere (upper)
    T_atm =228.65+((h -32e3)*2.8e-3);
    P_atm =868*(T_atm /228.65)^(g/(R*-2.8e-3));
elseif h >=47e3 && h <=51e3 %stratopause
    T_atm =270.65;
     P_atm =111*exp((-1*g*(h -47e3)/(R*T_atm )));
elseif h >51e3 && h <71e3 %mesosphere (lower)
    T_atm =270.65-((h -51e3)*2.8e-3);
    P_atm =67*(T_atm /270.65)^(g/(R*(2.8e-3)));
elseif h >= 71e3 && h <84e3 %mesosphere (upper)
    T_atm =214.65-((h -71e3)*2e-3);
    P_atm =3*(T_atm /214.65)^(g/(R*(2e-3)));
else 
     T_atm =214.65-((84e3-71e3)*2e-3); %space(vacuum)
     P_atm =0;
end 

rho_atm=P_atm/(R*T_atm);
end

function data =sim(sol)

global g vehicle t_seperation

stage=vehicle.stage;

m_pld=sol(1);
theta_kick=sol(2);
t_pitchover=sol(3);
t_pitch=sol(4);
t_burn_2=sol(5);
theta_SECO=sol(6);

opts_ode = odeset();
   
% Burn 1a (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1= -stage(1).m_prop / m_dot_1;
t_burn_1a=t_pitchover+t_pitch;

Y0_1a=[0;0;m0;0;pi/2]; %h=0, v=0, m=m0, x=0, theta=90degs

[t1a,Y1a]=ode45(@(t,Y1a) burn_1a(t,Y1a,stage,m_dot_1,g,t_pitchover,theta_kick,t_pitch),[0 t_burn_1a],Y0_1a,opts_ode);

% Burn 1b (gravity turn)

t_burn_1b=t_burn_1-t_burn_1a;
Y0_1b=Y1a(end,:)'; %h=0, v=0, m=m0, x=0, theta=90

[t1b,Y1b]=ode45(@(t,Y1b) burn_1b(t,Y1b,stage,m_dot_1,g),[t_burn_1a t_burn_1b+t_burn_1a],Y0_1b,opts_ode);

h_1b=Y1b(:,1);
v_1b=Y1b(:,2);
m_1b=Y1b(:,3);
x_1b=Y1b(:,4);
theta_1b=Y1b(:,5);

h_MECO = h_1b(end);
v_MECO = v_1b(end);
m_MECO = m_1b(end) - stage(1).m_dry;
x_MECO = x_1b(end);
theta_MECO = theta_1b(end);

% stage seperation

t_ignition_2=t_burn_1+t_seperation;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

h_sep=Y_sep(:,1);
h_SES=h_sep(end);

v_sep=Y_sep(:,2);
v_SES=v_sep(end);

m_sep=Y_sep(:,3);
m_SES=m_sep(end);

x_sep=Y_sep(:,4);
x_SES=x_sep(end);

theta_sep=Y_sep(:,5);
theta_SES=theta_sep(end);


% Burn 2 (up to apogee) 

m_dot_2= -stage(2).thrust/(stage(2).isp*g);

t_SECO_a=t_burn_2+t_burn_1+t_seperation;

Y0_2=[h_SES;v_SES;m_SES;x_SES;theta_SES]; %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,g),[t_ignition_2 t_SECO_a],Y0_2,opts_ode);


data.t=    [t1a      ;t1b      ;t_sep      ;t2     ];
data.h=    [Y1a(:,1) ;Y1b(:,1) ;Y_sep(:,1) ;Y2(:,1)];
data.v=    [Y1a(:,2) ;Y1b(:,2) ;Y_sep(:,2) ;Y2(:,2)];
data.m=    [Y1a(:,3) ;Y1b(:,3) ;Y_sep(:,3) ;Y2(:,3)];
data.x=    [Y1a(:,4) ;Y1b(:,4) ;Y_sep(:,4) ;Y2(:,4)];
data.theta=[Y1a(:,5) ;Y1b(:,5) ;Y_sep(:,5) ;Y2(:,5)];

end 

function [h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, ~ ,sol)

 global mu r_earth v_earth g vehicle

stage=vehicle.stage;
m_pld=sol(1);

h_SECO=data.h(end);
v_SECO=data.v(end);
theta_SECO=data.theta(end);
m_SECO=data.m(end);

r0 = [r_earth + h_SECO; 0];
v0 = [v_SECO*sin(theta_SECO); v_SECO*cos(theta_SECO) + v_earth]; % inertial

r0mag = norm(r0);
v0mag = norm(v0);

eps  = v0mag^2/2 - mu/r0mag;
hvec = r0(1)*v0(2) - r0(2)*v0(1);
habs = abs(hvec);

e = sqrt(1 + 2*eps*habs^2/mu^2);
a = -mu/(2*eps);

r_per = a*(1 - e);
r_apo = a*(1 + e);

h_per = r_per - r_earth;
h_apo = r_apo - r_earth;

% circularization at apogee (optional)
va = sqrt(mu*(2/r_apo - 1/a));
vcirc = sqrt(mu/r_apo);
delta_v_2b = max(vcirc - va, 0);

m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage(2).isp*g)));
m_prop_remaining=m_SECO-m_pld-stage(2).m_dry;

end


%% ODE functions

function dY1adt= burn_1a(t,Y1a,stage,m_dot_1,g,t_pitchover,theta_kick,t_pitch)
h=Y1a(1);
v=Y1a(2);
m=Y1a(3);
x=Y1a(4);
theta=Y1a(5);

%m=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_1*t;
Cd=stage(1).Cd;
A=stage(1).A;

[rho_atm,P_atm]=isa_model(h);

T=stage(1).thrust-(P_atm*(stage(1).A_exit));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;

if t<t_pitchover
dthetadt=0;
else
dthetadt= -(theta_kick/t_pitch) ;
end

dY1adt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

function dY1bdt= burn_1b(~,Y1b,stage,m_dot_1,g)

h=Y1b(1);
v=Y1b(2);
m=Y1b(3);
x=Y1b(4);
theta=Y1b(5);

Cd=stage(1).Cd;
A=stage(1).A;

[rho_atm,P_atm]=isa_model(h);

T=stage(1).thrust-(P_atm*(stage(1).A_exit));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY1bdt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

function dY_sepdt= seperation(~,Y_sep,stage,g)
h=Y_sep(1);
v=Y_sep(2);
m=Y_sep(3);
x=Y_sep(4);
theta=Y_sep(5);

T=0;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=0;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY_sepdt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

function dY2dt= burn_2(t,Y2,stage,m_dot_2,g)
h=Y2(1);
v=Y2(2);
m=Y2(3);
x=Y2(4);
theta=Y2(5);

%t_rel=t-t_burn_1;
T=stage(2).thrust;
%m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);


dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_2;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);


dY2dt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

    end
