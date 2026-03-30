
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

h_orbit=400e3;
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

azimuth=deg2rad(90); %due east 

vehicle=Electron();
stage=vehicle.stage;

m_pld=200;
t_coast=5;
theta_kick=deg2rad(4);
t_pitchover=10;
t_pitch=15;
t_burn_2a=350;


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
m_MECO=m_1b(end);
x_MECO=x_1b(end);
theta_MECO=theta_1b(end);

%% concatenate burn 1

t1=[t1a; t1b];
h_1=[h_1a; h_1b];
v_1=[v_1a; v_1b];
m_1=[m_1a; m_1b];
x_1=[x_1a; x_1b];
theta_1=[theta_1a; theta_1b];

%% Optimiser stage seperation

t_ignition_2a=t_burn_1+t_coast;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g),[t_burn_1 t_ignition_2a],Y0_sep);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
x_sep=Y_sep(:,4);
theta_sep=Y_sep(:,5);

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO-stage(1).m_dry;
x_SES=x_sep(end);
theta_SES=theta_sep(end);

%% Optimiser 2nd stage burn 

m_dot_2a= -stage(2).thrust/(stage(2).isp*g);

t_SECO=t_burn_2a+t_burn_1+t_coast;

Y0_2a=[h_SES;v_SES;m_SES;x_SES;theta_SES]; %state at SES

[t_2a,Y2a]=ode45(@(t,Y2a) burn_2(t,Y2a,stage,m_dot_2a,g),[t_ignition_2a t_SECO],Y0_2a);

h_2a=Y2a(:,1);
v_2a=Y2a(:,2);
m_2a=Y2a(:,3);
x_2a=Y2a(:,4);
theta_2a=Y2a(:,5);

v_SECO=v_2a(end);
h_SECO=h_2a(end);
x_SECO=x_2a(end);
m_SECO=m_2a(end);
theta_SECO=theta_2a(end);


%% Trajectory ODEs

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
T = max(T,0);

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

T=stage(1).thrust-(P_atm*(stage(1).A_exit));
T = max(T,0);

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

function dY2adt= burn_2(~,Y2a,stage,m_dot_2a,g)
h=Y2a(1);
v=Y2a(2);
m=Y2a(3);
x=Y2a(4);
theta=Y2a(5);

%t_rel=t-t_burn_1;
T=stage(2).thrust;
%m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);


dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_2a;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY2adt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end



    %% Dynamic pressure 

rho_1=zeros(size(h_1));
for i=1:length(h_1)
    h=h_1(i);

    [rho_1(i),~]=isa_model(h);
end

rho_2=zeros(size(h_2a));
for i=1:length(h_2a)
    h=h_2a(i);

    [rho_2(i),~]=isa_model(h);
end

q_1=(0.5.*rho_1.*v_1.^2)./1000; %units in kpa
q_2=(0.5.*rho_2.*v_2a.^2)./1000; %units in kpa
max_q=max(q_1);

[rho_atm,P_atm]=isa_model(h);

%% Acceleration

a_1=(stage(1).thrust-(1000.*q_1.*stage(1).A*stage(1).Cd)-m_1.*g)./(m_1.*g); %units in g
a_2=(stage(2).thrust-(1000.*q_2.*stage(2).A*stage(2).Cd)-m_2a.*g)./(m_2a.*g); %units in g

thrust_2=zeros(size(t_2a));
thrust_2=thrust_2+stage(2).thrust;
drag_2=1000.*q_2.*stage(2).A*stage(2).Cd;
weight_2=m_2a.*g;

%% Repeat orbital calcs

v_tan_inertial = v_SECO*cos(theta_SECO) + v_earth;
v_rad_inertial=v_SECO*sin(theta_SECO);
v_inertial=sqrt(v_rad_inertial^2+v_tan_inertial^2);

orbital_energy=(v_inertial^2)/2 - mu/(r_earth+h_SECO);
ang_mom=(r_earth+h_SECO)*v_tan_inertial;
e=sqrt(1+(2*orbital_energy*ang_mom^2)/mu^2)
a_semimajor=-mu/(2*orbital_energy);

r_apogee= a_semimajor*(1+e);
r_perigee=a_semimajor*(1-e);

h_apogee=r_apogee-r_earth
h_perigee=r_perigee-r_earth
h_min=120e3;

v_apogee=ang_mom/r_apogee
v_circ=sqrt(mu/r_apogee);

delta_v_2b=max(v_circ - v_apogee, 0);
m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage(2).isp*g)))
m_prop_remaining=m_SECO-m_pld-stage(2).m_dry


%% plots 

%converting to km
h_orbit_km=h_orbit./1000;
h_1_km=h_1./1000;
h_sep_km=h_sep./1000;
h_2_km=h_2a./1000;
x_1_km=x_1./1000;
x_sep_km=x_sep./1000;
x_2_km=x_2a./1000;
%converting to degrees
theta_1_deg=rad2deg(theta_1);
theta_sep_deg=rad2deg(theta_sep);
theta_2_deg=rad2deg(theta_2a);
theta_kick_deg=rad2deg(theta_kick);

figure;

subplot(2,2,1);
plot(t1,h_1_km) %altitude vs time
grid
hold on 
plot(t_sep,h_sep_km) 
plot(t_2a,h_2_km) 
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
legend("Burn 1",Location="northeast")
hold off 

subplot(2,2,3);
plot(x_1_km,h_1_km) %altitude vs downrange 
grid
hold on
plot(x_sep_km,h_sep_km)
plot(x_2_km,h_2_km)
xlabel("downrange distance (km)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off 

subplot(2,2,4);
plot(t1,theta_1_deg) % theta vs time
grid
hold on
plot (t_sep,theta_sep_deg)
plot (t_2a,theta_2_deg)
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

fprintf('t_burn_2 = %.2f s\n', t_burn_2a)
fprintf('theta_SECO = %.2f degrees\n', rad2deg(theta_SECO))
fprintf('max_q = %.2f kPa\n', max_q)

fprintf('t_coast = %.2f s\n', t_coast)
fprintf('theta_kick = %.2f degs\n', theta_kick_deg)
fprintf('t_pitchover = %.2f s\n', t_pitchover)
fprintf('t_pitch = %.2f s\n', t_pitch)
fprintf('m_pld = %.2f kg\n', m_pld)

