
% Proof of concept script 
clc
clear

%% Setup

R=287; %units=J/kg*K
g=9.80665; %units=m/s^2
m_pld=20000;

stage=[
    struct('m_dry', 25000, 'm_prop', 400000, 'isp', 283, 'thrust', 7.6e6, ...
    'Cd', 0.3, 'A', 10)
    struct('m_dry', 3900, 'm_prop', 92000, 'isp', 348, 'thrust', 981e3, ...
    'Cd', 0.25, 'A', 10) 
    ];

%% Burn 1

m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

Y0_1=[0;0]; %h=0, v=0

[t1,Y1]=ode45(@(t,Y1) burn_1(t,Y1,stage,m_dot_1,m_pld,R,g),[0 t_burn_1],Y0_1);

h_1=Y1(:,1);
v_1=Y1(:,2);

v_MECO=v_1(end);
h_MECO=h_1(end);

function dY1dt= burn_1(t,Y1,stage,m_dot_1,m_pld,R,g)
h=Y1(1);
v=Y1(2);

T=stage(1).thrust;
m=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_1*t;
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

dhdt=v;
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g)/m;

dY1dt = [dhdt;dvdt];

end

%% 2nd stage burn 

m_dot_2=-stage(2).thrust/(stage(2).isp * g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_end=t_burn_2+t_burn_1;

Y0_2=[h_MECO;v_MECO]; %h=height at seperation, v=velocity at seperation

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,m_pld,R,g,t_burn_1),[t_burn_1 t_end],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);

v_SECO=v_2(end);
h_SECO=h_2(end);

function dY2dt= burn_2(t,Y2,stage,m_dot_2,m_pld,R,g,t_burn_1)
h=Y2(1);
v=Y2(2);

t_rel=t-t_burn_1;
T=stage(2).thrust;
m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
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

dY2dt = [dhdt;dvdt];

end

%% plots 

h_1_km=h_1./1000;
h_2_km=h_2./1000;

plot(t1,h_1_km) 
grid
hold on 
plot(t2,h_2_km) 
xlabel("time (s)")
ylabel("altitude (km)")
legend("Burn 1","Burn 2")
hold off

figure
plot(t1,v_1)
grid
hold on
plot (t2,v_2)
xlabel("time (s)")
ylabel("velocity (m/s)")
legend("Burn 1","Burn 2")
hold off 


%% print outputs 

fprintf('h_MECO = %.2f km\n', h_MECO/1000)
fprintf('v_MECO = %.2f m/s\n', v_MECO)
fprintf('t_burn_1 = %.2f s\n', t_burn_1)
fprintf('h_SECO = %.2f km\n', h_SECO/1000)
fprintf('v_SECO = %.2f m/s\n', v_SECO)
fprintf('t_burn_2 = %.2f s\n', t_burn_2)