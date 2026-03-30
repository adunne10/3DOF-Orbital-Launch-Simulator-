
close all
clc
clear

%% Setup

g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m

h_orbit=220e3;
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));
theta_SECO_goal=deg2rad(0.1);

stage=[
    struct('m_dry', 25000, 'm_prop', 400000, 'isp', 283, 'thrust', 7.77e6, ...
    'Cd', 0.3, 'A', 10, 'A_exit', 0.9)
    struct('m_dry', 3900, 'm_prop', 92670, 'isp', 348, 'thrust', 981e3, ...
    'Cd', 0.25, 'A', 10, 'A_exit',0) 
    ];

%% Optimisation

sol_0=   [10000 ; 5  ; deg2rad(3) ; 10; 7]; %initial guess for m_pld , t_coast , theta_kick , t_pitchover & t_pitch
lowbound=[ 5000 ; 4 ; deg2rad(2) ; 5 ; 5]; %best guess for lower bounds
upbound= [ 22000 ; 6 ; deg2rad(10) ; 20 ; 10]; %best guess for upper bounds

objective=@(sol) -sol(1)/1e3;
nonlcon=@(sol) constraints(sol,stage,g,h_orbit,v_orbit,m_earth,r_earth,theta_SECO_goal);

options=optimoptions('fmincon', 'Algorithm', 'sqp', 'MaxFunctionEvaluations',1000 ,'Display','iter','StepTolerance',1e-8,'ConstraintTolerance',1e-4, 'MaxIterations', 1000);

problem=createOptimProblem('fmincon','x0',sol_0,'objective',objective,'lb',lowbound,'ub',upbound,'nonlcon',nonlcon,'options',options);
ms=MultiStart('UseParallel',false,'Display','iter');

num_starts=10;
rand_starts=rand(num_starts,length(sol_0)).*(upbound-lowbound)'+lowbound';
start_points=CustomStartPointSet(rand_starts);

[sol_optimal,fval,exitflag,output,solutions]=run(ms,problem,start_points);

m_pld_max=sol_optimal(1);
t_coast_optimal=sol_optimal(2);
theta_kick_optimal=sol_optimal(3);
t_pitchover_optimal=sol_optimal(4);
t_pitch_optimal=sol_optimal(5);

function [c,ceq] = constraints(sol,stage,g,h_target,v_target,~,~,theta_SECO_goal)

m_pld=sol(1);
t_coast=sol(2);
theta_kick=sol(3);
t_pitchover=sol(4);
t_pitch=sol(5);

[h_SECO, v_SECO] = trajectory(m_pld,stage,g,t_coast,t_pitchover, theta_kick,theta_SECO_goal,t_pitch);

    % equality constraints: hit target altitude, velocity, horizontal theta=0
    ceq =  [
            (h_SECO - h_target)/1000;
            (v_SECO - v_target)/100 
           ];
    c = [];
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

function [h_SECO, v_SECO] = trajectory( m_pld,stage,g,t_coast,t_pitchover, theta_kick,theta_SECO_goal,t_pitch)

%% Optimiser Burn 1a (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

t_burn_1a=t_pitchover+t_pitch;

Y0_1a=[0;0;m0;0;pi/2]; %h=0, v=0, m=m0, x=0, theta=90

[t1a,Y1a]=ode45(@(t,Y1a) burn_1a(t,Y1a,stage,m_dot_1,g,t_pitchover,theta_kick,t_pitch),[0 t_burn_1a],Y0_1a);

h_1a=Y1a(:,1);
v_1a=Y1a(:,2);
m_1a=Y1a(:,3);
x_1a=Y1a(:,4);
theta_1a=Y1a(:,5);

v_pitchover=v_1a(end);
h_pitchover=h_1a(end);
m_pitchover=m_1a(end);
x_pitchover=x_1a(end);
theta_pitchover=theta_1a(end);

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

T=stage(1).thrust+(P_atm*(stage(1).A_exit));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;

if t<=t_pitchover
dthetadt=0;
else
dthetadt=-theta_kick/t_pitch;
end

dY1adt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% Optimiser Burn 1b (gravity turn)

t_burn_1b=t_burn_1-t_burn_1a;

Y0_1b=[h_pitchover;v_pitchover;m_pitchover;x_pitchover;theta_pitchover]; %h=0, v=0, m=m0, x=0, theta=90

[t1b,Y1b]=ode45(@(t,Y1b) burn_1b(t,Y1b,stage,m_dot_1,g),[t_burn_1a t_burn_1b+t_burn_1a],Y0_1b);

h_1b=Y1b(:,1);
v_1b=Y1b(:,2);
m_1b=Y1b(:,3);
x_1b=Y1b(:,4);
theta_1b=Y1b(:,5);

v_MECO=v_1b(end);
h_MECO=h_1b(end);
m_MECO=stage(2).m_dry + stage(2).m_prop + m_pld;
x_MECO=x_1b(end);
theta_MECO=theta_1b(end);

function dY1bdt= burn_1b(~,Y1b,stage,m_dot_1,g)

h=Y1b(1);
v=Y1b(2);
m=Y1b(3);
x=Y1b(4);
theta=Y1b(5);

%m=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_1*t;
Cd=stage(1).Cd;
A=stage(1).A;

[rho_atm,P_atm]=isa_model(h);

T=stage(1).thrust+(P_atm*(stage(1).A_exit));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY1bdt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% Optimiser stage seperation

t_ignition_2=t_burn_1+t_coast;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[~,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g),[t_burn_1 t_ignition_2],Y0_sep);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
x_sep=Y_sep(:,4);
theta_sep=Y_sep(:,5);

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO;
x_SES=x_sep(end);
theta_SES=theta_sep(end);

    function dY_sepdt= seperation(~,Y_sep,stage,g)
h=Y_sep(1);
v=Y_sep(2);
m=Y_sep(3);
x=Y_sep(4);
theta=Y_sep(5);

T=0;
%m=stage(2).m_dry + stage(2).m_prop + m_pld;
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

%% Optimiser 2nd stage burn 

m_dot_2= -stage(2).thrust/(stage(2).isp*g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_SECO=t_burn_2+t_burn_1+t_coast;

dthetadt_SES=(-(g*((6371e3/(6371e3+h_SES))^2))/(v_SES+0.1)-((v_SES+0.1)/(h_SES+6371e3)))*cos(theta_SES);

a = [1 0 0 0;
     1 1 1 1;
     0 1 0 0;
     0 1 2 3];

b = [theta_SES; theta_SECO_goal; dthetadt_SES*t_burn_2; 0];

theta_coeff = a\b;

Y0_2=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,g,t_ignition_2,t_burn_2,theta_coeff),[t_ignition_2 t_SECO],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);
x_2=Y2(:,4);

    tau_2 = (t2 - t_ignition_2)./t_burn_2;
    tau_2 = min(max(tau_2,0),1);

    % Evaluate cubic theta and derivative
    theta_2 = theta_coeff(1) + theta_coeff(2).*tau_2 + theta_coeff(3).*tau_2.^2 + theta_coeff(4).*tau_2.^3;

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
theta_SECO=theta_2(end);

    function dY2dt= burn_2(t,Y2,stage,m_dot_2,g,t_ignition_2 ,t_burn_2,theta_coeff)
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

    tau = (t - t_ignition_2)/t_burn_2;
    tau = min(max(tau,0),1);

    % Evaluate cubic theta and derivative
    theta = theta_coeff(1) + theta_coeff(2)*tau + theta_coeff(3)*tau.^2 + theta_coeff(4)*tau.^3;
    dthetadt = (theta_coeff(2) + 2*theta_coeff(3)*tau + 3*theta_coeff(4)*tau.^2)/t_burn_2;

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_2;

dY2dt = [dhdt;dvdt;dmdt;dxdt];

end


end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Burn 1a (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld_max;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

t_burn_1a=t_pitchover_optimal+t_pitch_optimal;

Y0_1a=[0;0;m0;0;pi/2]; %h=0, v=0, m=m0, x=0, theta=90

[t1a,Y1a]=ode45(@(t,Y1a) burn_1a(t,Y1a,stage,m_dot_1,g,t_pitchover_optimal,theta_kick_optimal,t_pitch_optimal),[0 t_burn_1a],Y0_1a);

h_1a=Y1a(:,1);
v_1a=Y1a(:,2);
m_1a=Y1a(:,3);
x_1a=Y1a(:,4);
theta_1a=Y1a(:,5);

v_pitchover=v_1a(end);
h_pitchover=h_1a(end);
m_pitchover=m_1a(end);
x_pitchover=x_1a(end);
theta_pitchover=theta_1a(end);

    function dY1adt= burn_1a(t,Y1a,stage,m_dot_1,g,t_pitchover_optimal,theta_kick_optimal,t_pitch_optimal)
h=Y1a(1);
v=Y1a(2);
m=Y1a(3);
x=Y1a(4);
theta=Y1a(5);

%m=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_1*t;
Cd=stage(1).Cd;
A=stage(1).A;

[rho_atm,P_atm]=isa_model(h);

T=stage(1).thrust+(P_atm*(stage(1).A_exit));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;

if t<t_pitchover_optimal
dthetadt=0;
else
dthetadt= -(theta_kick_optimal/t_pitch_optimal) ;
end

dY1adt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% Burn 1b (gravity turn)

t_burn_1b=t_burn_1-t_burn_1a;

Y0_1b=[h_pitchover;v_pitchover;m_pitchover;x_pitchover;theta_pitchover]; %h=0, v=0, m=m0, x=0, theta=90

[t1b,Y1b]=ode45(@(t,Y1b) burn_1b(t,Y1b,stage,m_dot_1,g),[t_burn_1a t_burn_1b+t_burn_1a],Y0_1b);

h_1b=Y1b(:,1);
v_1b=Y1b(:,2);
m_1b=Y1b(:,3);
x_1b=Y1b(:,4);
theta_1b=Y1b(:,5);

v_MECO=v_1b(end);
h_MECO=h_1b(end);
m_MECO=stage(2).m_dry + stage(2).m_prop + m_pld_max;
x_MECO=x_1b(end);
theta_MECO=theta_1b(end);

function dY1bdt= burn_1b(~,Y1b,stage,m_dot_1,g)

h=Y1b(1);
v=Y1b(2);
m=Y1b(3);
x=Y1b(4);
theta=Y1b(5);

Cd=stage(1).Cd;
A=stage(1).A;

[rho_atm,P_atm]=isa_model(h);

T=stage(1).thrust+(P_atm*(stage(1).A_exit));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY1bdt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end


%% concatenate burn 1

t1=[t1a; t1b];
h_1=[h_1a; h_1b];
v_1=[v_1a; v_1b];
m_1=[m_1a; m_1b];
x_1=[x_1a; x_1b];
theta_1=[theta_1a; theta_1b];

%% stage seperation

t_ignition_2=t_burn_1+t_coast_optimal;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g),[t_burn_1 t_ignition_2],Y0_sep);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
m_sep=Y_sep(:,3);
x_sep=Y_sep(:,4);
theta_sep=Y_sep(:,5);

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO;
x_SES=x_sep(end);
theta_SES=theta_sep(end);

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

%% 2nd stage burn 

m_dot_2=-stage(2).thrust/(stage(2).isp*g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_SECO=t_burn_2+t_burn_1+t_coast_optimal;

dthetadt_SES=(-(g*((6371e3/(6371e3+h_SES))^2))/(v_SES+0.1)-((v_SES+0.1)/(h_SES+6371e3)))*cos(theta_SES);

a = [1 0 0 0;
     1 1 1 1;
     0 1 0 0;
     0 1 2 3];

b = [theta_SES; theta_SECO_goal; dthetadt_SES*t_burn_2; 0];

theta_coeff = a\b;

   
Y0_2=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,g,t_ignition_2,t_burn_2,theta_coeff),[t_ignition_2 t_SECO],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);
x_2=Y2(:,4);

    tau_2 = (t2 - t_ignition_2)./t_burn_2;
    tau_2 = min(max(tau_2,0),1);

    % Evaluate cubic theta and derivative
    theta_2 = theta_coeff(1) + theta_coeff(2).*tau_2 + theta_coeff(3).*tau_2.^2 + theta_coeff(4).*tau_2.^3;

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
theta_SECO=theta_2(end);

function dY2dt= burn_2(t,Y2,stage,m_dot_2,g,t_ignition_2,t_burn_2,theta_coeff)
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

    tau = (t - t_ignition_2)/t_burn_2;
    tau = min(max(tau,0),1);

    % Evaluate cubic theta and derivative
    theta = theta_coeff(1) + theta_coeff(2)*tau + theta_coeff(3)*tau.^2 + theta_coeff(4)*tau.^3;
    dthetadt = (theta_coeff(2) + 2*theta_coeff(3)*tau + 3*theta_coeff(4)*tau.^2)/t_burn_2;

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

a_1=(stage(1).thrust-(1000.*q_1.*stage(1).A*stage(1).Cd)-m_1.*g)./(m_1.*g); %units in g
a_2=(stage(2).thrust-(1000.*q_2.*stage(2).A*stage(2).Cd)-m_2.*g)./(m_2.*g); %units in g

thrust_2=zeros(size(t2));
thrust_2=thrust_2+stage(2).thrust;
drag_2=1000.*q_2.*stage(2).A*stage(2).Cd;
weight_2=m_2.*g;

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
theta_1_deg=rad2deg(abs(theta_1-pi/2));
theta_sep_deg=rad2deg(abs(theta_sep-pi/2));
theta_2_deg=rad2deg(abs(theta_2-pi/2));
theta_kick_optimal_deg=rad2deg(theta_kick_optimal);

figure;

subplot(2,2,1);
plot(t1,h_1_km) %altitude vs time
grid
hold on 
plot(t_sep,h_sep_km) 
plot(t2,h_2_km) 
xlabel("time (s)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off

subplot(2,2,2);
plot(x_1_km,h_1_km) %altitude vs downrange distance
grid
hold on 
plot(x_sep_km,h_sep_km) 
plot(x_2_km,h_2_km) 
xlabel("downrange distance (km)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off

subplot(2,2,3);
plot(t1,v_1) % veolcity vs time
grid
hold on
plot (t_sep,v_sep)
plot (t2,v_2)
ylabel("velocity (m/s)")
xlabel("time (s)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off 

subplot(2,2,4);
plot(t1,theta_1_deg) % theta vs time
grid
hold on
plot (t_sep,theta_sep_deg)
plot (t2,theta_2_deg)
xlabel("time (s)")
ylabel("theta (deg)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off 


%% print outputs 
 
fprintf('h_MECO = %.2f km\n', h_MECO/1000)
fprintf('v_MECO = %.2f m/s\n', v_MECO)
fprintf('t_burn_1 = %.2f s\n', t_burn_1)
fprintf('h_SECO = %.2f km\n', h_SECO/1000)
fprintf('v_SECO = %.2f m/s\n', v_SECO)
fprintf('theta_SECO = %.2f degrees\n', rad2deg(theta_SECO))
fprintf('t_burn_2 = %.2f s\n', t_burn_2)
fprintf('max_q = %.2f kPa\n', max_q)
fprintf('v_orbit = %.2f m/s\n', v_orbit)
fprintf('h_orbit = %.2f km\n', h_orbit_km)

fprintf('t_coast_optimal = %.2f s\n', t_coast_optimal)
fprintf('theta_kick_optimal = %.2f degs\n', theta_kick_optimal_deg)
fprintf('t_pitchover_optimal = %.2f s\n', t_pitchover_optimal)
fprintf('t_pitch_optimal = %.2f s\n', t_pitch_optimal)
fprintf('m_pld_max = %.2f kg\n', m_pld_max)

