
close all
clc
clear

%% Setup

R=287; %units=J/kg*K
g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m

h_orbit=350e3;
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

stage=[
    struct('m_dry', 25000, 'm_prop', 400000, 'isp', 283, 'thrust', 7.6e6, ...
    'Cd', 0.3, 'A', 10, 'A_exit', 0.9)
    struct('m_dry', 3900, 'm_prop', 92000, 'isp', 348, 'thrust', 981e3, ...
    'Cd', 0.25, 'A', 10, 'A_exit',0) 
    ];

alpha=-4.1764;
beta=4.9990;
charlie=-2.3672;
m_pld=5820;
t_coast=5;

%% Burn 1

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

Y0_1=[0;0;m0;0]; %h=0, v=0, m=m0, x=0, theta=pi/2 (vertical)

[t1,Y1]=ode45(@(t,Y1) burn_1(t,Y1,stage,m_dot_1,R,g,alpha,beta,charlie,h_orbit),[0 t_burn_1],Y0_1);

h_1=Y1(:,1);
v_1=Y1(:,2);
m_1=Y1(:,3);
x_1=Y1(:,4);
theta_1=pi/2 + alpha.*(h_1./h_orbit) + beta.*(h_1./h_orbit).^2 + charlie.*(h_1./h_orbit).^3;

v_MECO=v_1(end);
h_MECO=h_1(end);
m_MECO=stage(2).m_dry + stage(2).m_prop + m_pld;
x_MECO=x_1(end);
theta_MECO=theta_1(end);

    function dY1dt= burn_1(t,Y1,stage,m_dot_1,R,g,alpha,beta,charlie,h_orbit)
h=Y1(1);
v=Y1(2);
m=Y1(3);
x=Y1(4);


%m=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_1*t;
Cd=stage(1).Cd;
A=stage(1).A;

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

T=stage(1).thrust+(P_atm*(stage(1).A_exit));

theta=pi/2+alpha*(h/h_orbit)+beta*(h/h_orbit)^2+charlie*(h/h_orbit)^3;

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;

dY1dt = [dhdt;dvdt;dmdt;dxdt];

end

%% stage seperation

t_ignition_2=t_burn_1+t_coast;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,R,g,alpha,beta,charlie,h_orbit),[t_burn_1 t_ignition_2],Y0_sep);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
x_sep=Y_sep(:,4);
theta_sep=pi/2 + alpha.*(h_sep./h_orbit) + beta.*(h_sep./h_orbit).^2 + charlie.*(h_sep./h_orbit).^3;

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO;
x_SES=x_sep(end);
theta_SES=theta_sep(end);

    function dY_sepdt= seperation(t,Y_sep,stage,R,g,alpha,beta,charlie,h_orbit)
h=Y_sep(1);
v=Y_sep(2);
m=Y_sep(3);
x=Y_sep(4);

T=0;
%m=stage(2).m_dry + stage(2).m_prop + m_pld;
Cd=stage(2).Cd;
A=stage(2).A;

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
     T_atm =214.65-((84e3-71e3)*2e-3); %space (vacuum)
     P_atm =0;
end 

rho_atm=P_atm/(R*T_atm);

theta=pi/2+alpha*(h/h_orbit)+beta*(h/h_orbit)^2+charlie*(h/h_orbit)^3;

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2))/m;
dmdt=0;

dY_sepdt = [dhdt;dvdt;dmdt;dxdt];

end

%% 2nd stage burn 

m_dot_2=-stage(2).thrust/(stage(2).isp*g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_SECO=t_burn_2+t_burn_1+t_coast;

Y0_2=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,R,g,alpha,beta,charlie,h_orbit),[t_ignition_2 t_SECO],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);
x_2=Y2(:,4);
theta_2=pi/2 +alpha.*(h_2./h_orbit)+beta.*(h_2./h_orbit).^2+charlie.*(h_2./h_orbit).^3;

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
theta_SECO=theta_2(end);

    function dY2dt= burn_2(t,Y2,stage,m_dot_2,R,g,alpha,beta,charlie,h_orbit)
h=Y2(1);
v=Y2(2);
m=Y2(3);
x=Y2(4);

%t_rel=t-t_burn_1;
T=stage(2).thrust;
%m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
Cd=stage(2).Cd;
A=stage(2).A;


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
     T_atm =214.65-((84e3-71e3)*2e-3); %space (vacuum)
     P_atm =0;
end 

rho_atm=P_atm/(R*T_atm);

theta=pi/2+alpha*(h/h_orbit)+beta*(h/h_orbit)^2+charlie*(h/h_orbit)^3;

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2))/m;
dmdt=m_dot_2;

dY2dt = [dhdt;dvdt;dmdt;dxdt];

end

%% Dynamic pressure 

rho_1=zeros(size(h_1));
for i=1:length(h_1)
    h=h_1(i);

    [T_atm,P_atm,rho_1(i)]=ISA(h,R,g);
end

rho_2=zeros(size(h_2));
for i=1:length(h_2)
    h=h_2(i);

    [T_atm,P_atm,rho_2(i)]=ISA(h,R,g);
end

q_1=(0.5.*rho_1.*v_1.^2)./1000; %units in kpa
q_2=(0.5.*rho_2.*v_2.^2)./1000; %units in kpa
max_q=max(q_1);

function[T_atm,P_atm,rho_atm]=ISA(h,R,g)

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
     T_atm =214.65-((84e3-71e3)*2e-3); %space (vacuum)
     P_atm =0;
end 

rho_atm=P_atm/(R*T_atm);
end 

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

plot(t1,h_1_km) %altitude vs time
grid
hold on 
plot(t_sep,h_sep_km) 
plot(t2,h_2_km) 
xlabel("time (s)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off

figure%2
plot(x_1_km,h_1_km) %altitude vs downrange distance
grid
hold on 
plot(x_sep_km,h_sep_km) 
plot(x_2_km,h_2_km) 
xlabel("downrange distance (km)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off

figure%3
plot(t1,v_1) % veolcity vs time
grid
hold on
plot (t_sep,v_sep)
plot (t2,v_2)
ylabel("velocity (m/s)")
xlabel("time (s)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off 


figure%3
plot(t1,theta_1_deg) % theta vs time
grid
hold on
plot (t_sep,theta_sep_deg)
plot (t2,theta_2_deg)
xlabel("time (s)")
ylabel("theta (deg)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off 

figure 
plot(q_1,h_1_km) %altitude vs dynamic pressure
grid
hold on
xlabel("dynamic pressure (kPa)")
ylabel("altitude (km)")
legend("Burn 1",Location="southeast")
xlim([0 60]);
ylim([0 30]);
hold off 
% 
% figure 
% plot(a_1,h_1_km) %altitude vs acceleration
% grid
% hold on
% xlabel("acceleration (g)")
% ylabel("altitude (km)")
% legend("Burn 1",Location="southeast")
% hold off 

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