
close all
clc
clear

%% Setup

g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m
mu=G*m_earth;

omega_earth=7.2921159e-5; %units=rads/s
lat_launch=deg2rad(28.5); %cape canaveral
v_earth=omega_earth*r_earth*cos(lat_launch); %pure eastward velocity

h_orbit=500e3;
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

azimuth=deg2rad(90); %due east 

vehicle=Atlas_V_401();
stage=vehicle.stage;

t_seperation=15;

burn_time_max=stage(2).m_prop/(stage(2).thrust/(stage(2).isp*g));
burn_time_guess=300+((burn_time_max-300)/2);
pld_guess=stage(1).pld_guess+((stage(2).pld_guess-stage(1).pld_guess)/2);


%% Optimisation

sol_0=   [ pld_guess          ; 0  ; 1   ;  burn_time_guess ]; %initial guess for m_pld , t_seperation , theta_kick , t_pitchover ,  t_pitch
lowbound=[ stage(1).pld_guess ; 0  ; 0   ;  300             ]; %best guess for lower bounds
upbound= [ stage(2).pld_guess ; 0  ; 10  ;  burn_time_max   ]; %best guess for upper bounds

objective=@(sol) -sol(1)/1e3;
nonlcon=@(sol) constraints(sol,stage,g,mu,h_orbit,r_earth,v_earth,t_seperation);

options=optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display','iter','StepTolerance',1e-8,'ConstraintTolerance',1e-4, 'MaxIterations', 1000);

problem=createOptimProblem('fmincon','x0',sol_0,'objective',objective,'lb',lowbound,'ub',upbound,'nonlcon',nonlcon,'options',options);
ms=MultiStart('UseParallel',false,'Display','iter');

num_starts=15; 
rand_starts=rand(num_starts,length(sol_0)).*(upbound-lowbound)'+lowbound';
start_points=CustomStartPointSet(rand_starts);

[sol_optimal,fval,exitflag,output,solutions]=run(ms,problem,start_points);


m_pld_max=sol_optimal(1);
alpha_optimal=sol_optimal(2);
beta_optimal=sol_optimal(3);
t_burn_2_optimal=sol_optimal(4);

function [c,ceq] = constraints(sol,stage,g,mu,h_target,r_earth,v_earth,t_seperation)

m_pld=sol(1);
alpha=sol(2);
beta=sol(3);
t_burn_2=sol(4);

[h_SECO, v_SECO, m_SECO, theta_SECO] = trajectory(m_pld,stage,g,t_burn_2,t_seperation,alpha,beta);

%% Orbit calcs
r0 = [r_earth + h_SECO; 0];
v0 = [v_SECO*sin(theta_SECO); v_SECO*cos(theta_SECO) + v_earth]; % inertial

r0mag = norm(r0);
v0mag = norm(v0);

eps  = v0mag^2/2 - mu/r0mag;
hvec = r0(1)*v0(2) - r0(2)*v0(1);
habs = abs(hvec);

e = sqrt(1 + 2*eps*habs^2/mu^2);
a = -mu/(2*eps);

rp = a*(1 - e);
ra = a*(1 + e);

hp = rp - r_earth;
ha = ra - r_earth;

% circularization at apogee 
va = sqrt(mu*(2/ra - 1/a));
vcirc = sqrt(mu/ra);
delta_v_2b = max(vcirc - va, 0);

h_min=120e3;
m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage(2).isp*g)));
m_prop_remaining=m_SECO-m_pld-stage(2).m_dry;
h_cut_min = 150e3;         

    % equality constraints: hit target altitude, mass
    ceq =  [(ha - h_target)/1000;
           ];
    c = [(m_prop_2b-m_prop_remaining)/10;
         (h_min-hp)/1000;
         (h_cut_min - h_SECO)/1000;
         -theta_SECO
    ];
end 

function[rho_atm,P_atm] = isa_model(h)

R=287; %units=J/kg*K
g=9.80665; %units=m/s^2

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

function [h_SECO,v_SECO,m_SECO, theta_SECO] = trajectory( m_pld,stage,g,t_burn_2,t_seperation,alpha,beta)

%% Optimiser Burn 1 (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;


Y0_1=[0;0;m0;0]; %h=0, v=0, m=m0, x=0, theta=90

opts_ode = odeset();
[t1,Y1]=ode45(@(t,Y1) burn_1(t,Y1,stage,m_dot_1,g,alpha,beta),[0 t_burn_1],Y0_1,opts_ode);

h_1=Y1(:,1);
v_1=Y1(:,2);
m_1=Y1(:,3);
x_1=Y1(:,4);

theta_1 = atan( 1./(alpha + beta*t1) );

v_MECO=v_1(end);
h_MECO=h_1(end);
m_MECO=m_1(end)-stage(1).m_dry;
x_MECO=x_1(end);
theta_MECO=theta_1(end);


%% Optimiser stage seperation

t_ignition_2=t_burn_1+t_seperation;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g,alpha,beta,t_burn_1),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
x_sep=Y_sep(:,4);


theta_sep = atan( 1./(alpha + beta*t_sep) );

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO;
x_SES=x_sep(end);
theta_SES=theta_sep(end);

%% Optimiser 2nd stage burn 

m_dot_2= -stage(2).thrust/(stage(2).isp*g);

t_SECO=t_burn_2+t_burn_1+t_seperation;

Y0_2=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t_2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,g,alpha,beta,t_ignition_2),[t_ignition_2 t_SECO],Y0_2,opts_ode);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);
x_2=Y2(:,4);

theta_2 = atan( 1./(alpha + beta*t_2) );

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
m_SECO=m_2(end);
theta_SECO=theta_2(end);

%% Optimiser Trajectory Functions

function dY1dt= burn_1(t,Y1,stage,m_dot_1,g,alpha,beta)
h=Y1(1);
v=Y1(2);
m=Y1(3);
x=Y1(4);

%m=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_1*t;
Cd=stage(1).Cd;
A=stage(1).A;

[rho_atm,P_atm]=isa_model(h);

T=stage(1).thrust-(P_atm*(stage(1).A_exit));
T = max(T,0);

theta = atan( 1./(alpha + beta*t) );

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;

dY1dt = [dhdt;dvdt;dmdt;dxdt];

end

    function dY_sepdt= seperation(t,Y_sep,stage,g,alpha,beta,t_burn_1)
h=Y_sep(1);
v=Y_sep(2);
m=Y_sep(3);
x=Y_sep(4);

T=0;
%m=stage(2).m_dry + stage(2).m_prop + m_pld;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);

theta = atan( 1./(alpha + beta*t) );

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=0;

dY_sepdt = [dhdt;dvdt;dmdt;dxdt];

end

    function dY2dt= burn_2(t,Y2,stage,m_dot_2,g,alpha,beta,t_ignition_2)
h=Y2(1);
v=Y2(2);
m=Y2(3);
x=Y2(4);

%t_rel=t-t_burn_1;
T=stage(2).thrust;
%m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);

theta = atan( 1./(alpha + beta*t) );

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_2;

dY2dt = [dhdt;dvdt;dmdt;dxdt];

end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Burn 1 (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld_max;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1= -stage(1).m_prop / m_dot_1;

Y0_1=[0;0;m0;0]; %h=0, v=0, m=m0, x=0, theta=90

opts_ode = odeset();
[t1,Y1]=ode45(@(t,Y1) burn_1(t,Y1,stage,m_dot_1,g,alpha_optimal,beta_optimal),[0 t_burn_1],Y0_1,opts_ode);

h_1=Y1(:,1);
v_1=Y1(:,2);
m_1=Y1(:,3);
x_1=Y1(:,4);

theta_1 = atan( 1./(alpha_optimal + beta_optimal*t1) );

v_MECO=v_1(end);
h_MECO=h_1(end);
m_MECO=m_1(end);
x_MECO=x_1(end);
theta_MECO=theta_1(end);

%% stage seperation

t_ignition_2=t_burn_1+t_seperation;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g,alpha_optimal,beta_optimal,t_burn_1),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
m_sep=Y_sep(:,3);
x_sep=Y_sep(:,4);

theta_sep = atan( 1./(alpha_optimal + beta_optimal*t_sep) );

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO;
x_SES=x_sep(end);
theta_SES=theta_sep(end);
  

%% Burn 2 (up to apogee) 

m_dot_2= -stage(2).thrust/(stage(2).isp*g);

t_SECO_a=t_burn_2_optimal+t_burn_1+t_seperation;

dthetadt_SES=(-(g*((6371e3/(6371e3+h_SES))^2))/(v_SES+0.1)-((v_SES+0.1)/(h_SES+6371e3)))*cos(theta_SES);

Y0_2_a=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t_2,Y2_a]=ode45(@(t,Y2_a) burn_2(t,Y2_a,stage,m_dot_2,g,alpha_optimal,beta_optimal,t_ignition_2),[t_ignition_2 t_SECO_a],Y0_2_a,opts_ode);

h_2=Y2_a(:,1);
v_2=Y2_a(:,2);
m_2=Y2_a(:,3);
x_2=Y2_a(:,4);

theta_2 = atan( 1./(alpha_optimal + beta_optimal*t_2) );

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
m_SECO=m_2(end);
theta_SECO=theta_2(end);

%% redo orbit calcs
r0 = [r_earth + h_SECO; 0];
v0 = [v_SECO*sin(theta_SECO); v_SECO*cos(theta_SECO) + v_earth]; % inertial

r0mag = norm(r0);
v0mag = norm(v0);

eps  = v0mag^2/2 - mu/r0mag;
hvec = r0(1)*v0(2) - r0(2)*v0(1);
habs = abs(hvec);

e = sqrt(1 + 2*eps*habs^2/mu^2);
a = -mu/(2*eps);

rp = a*(1 - e);
ra = a*(1 + e);

hp = rp - r_earth;
ha = ra - r_earth;

% circularization at apogee (optional)
va = sqrt(mu*(2/ra - 1/a));
vcirc = sqrt(mu/ra);
delta_v_2b = max(vcirc - va, 0);

m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage(2).isp*g)));
m_prop_remaining=m_SECO-m_pld_max-stage(2).m_dry;

%% ODE functions

function dY1dt= burn_1(t,Y1,stage,m_dot_1,g,alpha_optimal,beta_optimal)
h=Y1(1);
v=Y1(2);
m=Y1(3);
x=Y1(4);

%m=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_1*t;
Cd=stage(1).Cd;
A=stage(1).A;

[rho_atm,P_atm]=isa_model(h);

T=stage(1).thrust-(P_atm*(stage(1).A_exit));

theta = atan( 1./(alpha_optimal + beta_optimal*t) );

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;


dY1dt = [dhdt;dvdt;dmdt;dxdt];

end


function dY_sepdt= seperation(t,Y_sep,stage,g,alpha_optimal,beta_optimal,t_burn_1)
h=Y_sep(1);
v=Y_sep(2);
m=Y_sep(3);
x=Y_sep(4);

T=0;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);

theta = atan( 1./(alpha_optimal + beta_optimal*t) );

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=0;

dY_sepdt = [dhdt;dvdt;dmdt;dxdt];

end

function dY2dt= burn_2(t,Y2,stage,m_dot_2,g,alpha_optimal,beta_optimal,t_ignition_2)
h=Y2(1);
v=Y2(2);
m=Y2(3);
x=Y2(4);

%t_rel=t-t_burn_1;
T=stage(2).thrust;
%m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);

theta = atan( 1./(alpha_optimal + beta_optimal*t) );

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_2;

dY2dt = [dhdt;dvdt;dmdt;dxdt];

    end

    %% Dynamic pressure 

rho_1=zeros(size(h_1));
for i=1:length(h_1)
    h=h_1(i);

    [rho_1(i),~]=isa_model(h);
end

rho_2=zeros(size(h_2));
for i=1:length(h_2)
    h=h_2(i);

    [rho_2(i),~]=isa_model(h);
end

q_1=(0.5.*rho_1.*v_1.^2)./1000; %units in kpa
q_2=(0.5.*rho_2.*v_2.^2)./1000; %units in kpa
max_q=max(q_1);

[rho_atm,P_atm]=isa_model(h);

%% Acceleration

g_local_1 = g*(r_earth/(r_earth+h_1)).^2;
g_local_2 = g*(r_earth/(r_earth+h_2)).^2;

a_1=(stage(1).thrust-(1000.*q_1.*stage(1).A*stage(1).Cd)-m_1.*g_local_1)./(m_1.*g_local_1); %units in g
a_2=(stage(2).thrust-(1000.*q_2.*stage(2).A*stage(2).Cd)-m_2.*g_local_2)./(m_2.*g_local_2); %units in g

thrust_2=zeros(size(t_2));
thrust_2=thrust_2+stage(2).thrust;
drag_2=1000.*q_2.*stage(2).A*stage(2).Cd;

%% plots 

%converting to km
h_orbit_km=h_orbit./1000;
h_1_km=h_1./1000;
h_sep_km=h_sep./1000;
h_2_km=h_2./1000;
x_1_km=x_1./1000;
x_sep_km=x_sep./1000;
x_2_km=x_2./1000;
%converting to degrees
theta_1_deg=rad2deg(theta_1);
theta_sep_deg=rad2deg(theta_sep);
theta_2_deg=rad2deg(theta_2);

figure;

subplot(2,2,1);
plot(t1,h_1_km) %altitude vs time
grid
hold on 
plot(t_sep,h_sep_km) 
plot(t_2,h_2_km) 
xlabel("time (s)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off

subplot(2,2,2);
plot(t1,q_1) %altitude vs dynamic pressure
grid
hold on
ylabel("dynamic pressure (kPa)")
xlabel("time (s)")
legend("Burn 1",Location="southeast")
hold off 

subplot(2,2,3);
plot(t1,a_1) %altitude vs acceleration
grid
hold on
xlabel("time (s)")
ylabel("acceleration (g)")
legend("Burn 1",Location="southeast")
hold off 

subplot(2,2,4);
plot(t1,theta_1_deg) % theta vs time
grid
hold on
plot (t_sep,theta_sep_deg)
plot (t_2,theta_2_deg)
xlabel("time (s)")
ylabel("theta (deg)")
legend("Burn 1","Stage seperation","Burn 2",Location="northeast")
hold off 


%% print outputs 
 
fprintf('h_MECO = %.2f km\n', h_MECO/1000)
fprintf('v_MECO = %.2f m/s\n', v_MECO)
fprintf('t_burn_1 = %.2f s\n', t_burn_1)
fprintf('h_SECO = %.2f km\n', h_SECO/1000)
fprintf('h_orbit = %.2f km\n', h_orbit_km)
fprintf('v_SECO = %.2f m/s\n', v_SECO)
fprintf('v_orbit = %.2f m/s\n', v_orbit)

fprintf('t_burn_2 = %.2f s\n', t_burn_2_optimal)
fprintf('theta_SECO = %.2f degrees\n', rad2deg(theta_SECO))
fprintf('max_q = %.2f kPa\n', max_q)


fprintf('m_prop_2b = %.2f kg\n', m_prop_2b)
fprintf('m_prop_remaining = %.2f kg\n', m_prop_remaining)

fprintf('h_apogee = %.2f km\n', ha/1000)
fprintf('h_perigee = %.2f km\n', hp/1000)
fprintf('e = %.2f \n', e)
fprintf('m_pld_max = %.2f kg\n', m_pld_max+(m_prop_remaining-m_prop_2b))
fprintf('alpha = %.2f \n', alpha_optimal)
fprintf('beta = %.2f \n', beta_optimal)