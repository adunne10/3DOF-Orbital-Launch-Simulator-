
% Proof of concept script 
close all
clc
clear


%% Setup

R=287; %units=J/kg*K
g=9.80665; %units=m/s^2
m_pld=20000;

stage=[
    struct('m_dry', 25000, 'm_prop', 400000, 'isp', 283, 'thrust', 7.6e6, ...
    'Cd', 0.3, 'A', 10, 'A_exit', 0.9)
    struct('m_dry', 3900, 'm_prop', 92000, 'isp', 348, 'thrust', 981e3, ...
    'Cd', 0.25, 'A', 10, 'A_exit',0) 
    ];

%% Burn 1

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

Y0_1=[0;0;m0]; %h=0, v=0

[t1,Y1]=ode45(@(t,Y1) burn_1(t,Y1,stage,m_dot_1,R,g),[0 t_burn_1],Y0_1);

h_1=Y1(:,1);
v_1=Y1(:,2);
m_1=Y1(:,3);

v_MECO=v_1(end);
h_MECO=h_1(end);

function dY1dt= burn_1(~,Y1,stage,m_dot_1,R,g)
h=Y1(1);
v=Y1(2);
m=Y1(3);

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
     T_atm =214.65-((84e3-71e3)*2e-3); %space (vacuum)
     P_atm =0;
end 

rho_atm=P_atm/(R*T_atm);
T=stage(1).thrust+(P_atm*(stage(1).A_exit));

dhdt=v;
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g)/m;
dmdt=m_dot_1;

dY1dt = [dhdt;dvdt;dmdt];

end
%% stage seperation

t_seperation=7; %varies from 5-10 seconds
t_ignition_2=t_burn_1+t_seperation;
m_MECO=stage(2).m_dry + stage(2).m_prop + m_pld;

Y0_sep=[h_MECO;v_MECO;m_MECO]; %h=height at seperation, v=velocity at seperation

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,R,g),[t_burn_1 t_ignition_2],Y0_sep);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);

v_ign=v_sep(end);
h_ign=h_sep(end);

function dY_sepdt= seperation(~,Y_sep,stage,R,g)
h=Y_sep(1);
v=Y_sep(2);
m=Y_sep(3);

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

dhdt=v;
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g)/m;
dmdt=0;

dY_sepdt = [dhdt;dvdt;dmdt];

end

%% 2nd stage burn 

m_dot_2=-stage(2).thrust/(stage(2).isp * g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_SECO=t_burn_2+t_burn_1+t_seperation;
m_ign=m_MECO; %no mass change 

Y0_2=[h_ign;v_ign;m_ign]; %h=height at ignition, v=velocity at ignition

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,R,g),[t_ignition_2 t_SECO],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);

v_SECO=v_2(end);
h_SECO=h_2(end);

function dY2dt= burn_2(~,Y2,stage,m_dot_2,R,g)
h=Y2(1);
v=Y2(2);
m=Y2(3);

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

dhdt=v;
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g)/m;
dmdt=m_dot_2;

dY2dt = [dhdt;dvdt;dmdt];

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

h_1_km=h_1./1000;
h_sep_km=h_sep./1000;
h_2_km=h_2./1000;

plot(t1,h_1_km) %altitude vs time
grid
hold on 
plot(t_sep,h_sep_km) 
plot(t2,h_2_km) 
xlabel("time (s)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2")
hold off

figure 
plot(v_1,h_1_km) % velocity vs altitude
grid
hold on
plot (v_sep,h_sep_km)
plot (v_2,h_2_km)
xlabel("velocity (m/s)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2")
hold off 

figure 
plot(v_1,h_1_km) % velocity vs altitidue
grid
hold on
plot (v_sep,h_sep_km)
plot (v_2,h_2_km)
xlabel("velocity (m/s)")
ylabel("altitude (km)")
legend("Burn 1","Stage seperation","Burn 2")
xlim([0 1500]);
ylim([0 30]);

hold off 

figure 
plot(q_1,h_1_km) %altitude vs dynamic pressure
grid
hold on
xlabel("dynamic pressure (kPa)")
ylabel("altitude (km)")
legend("Burn 1")
xlim([0 60]);
ylim([0 30]);
hold off 

figure 
plot(a_1,h_1_km) %altitude vs acceleration
grid
hold on

xlabel("acceleration (g)")
ylabel("altitude (km)")
legend("Burn 1")
hold off 


%% print outputs 

fprintf('h_MECO = %.2f km\n', h_MECO/1000)
fprintf('v_MECO = %.2f m/s\n', v_MECO)
fprintf('t_burn_1 = %.2f s\n', t_burn_1)
fprintf('h_SECO = %.2f km\n', h_SECO/1000)
fprintf('v_SECO = %.2f m/s\n', v_SECO)
fprintf('t_burn_2 = %.2f s\n', t_burn_2)
fprintf('max_q = %.2f kPa\n', max_q)