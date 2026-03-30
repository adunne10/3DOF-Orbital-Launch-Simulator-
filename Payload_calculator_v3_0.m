close all
clc
clear

%% Setup
%Constants
R=287; %units=J/kg*K
g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m

%target
h_orbit=400e3;
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

%LV characteristics
stage=[
    struct('m_dry', 25000, 'm_prop', 400000, 'isp', 283, 'thrust', 7.6e6, ...
    'Cd', 0.3, 'A', 10, 'A_exit', 0.9);
    struct('m_dry', 3900, 'm_prop', 92000, 'isp', 348, 'thrust', 981e3, ...
    'Cd', 0.25, 'A', 10, 'A_exit',0) 
    ];

sol_0=[deg2rad(2);500;10000]; %initial guess for theta_kick,h_pitchover and m_pld

lowbound=[deg2rad(0.1);50;10]; %best guess for lower bounds
upbound=[deg2rad(10);5000;50000]; %best guess for upper bounds

objective=@(sol) -sol(3);

nonlcon=@(sol) constraints(sol,stage,R,g,h_orbit,v_orbit,m_earth,r_earth);
options=optimoptions('fmincon', 'Display','iter','StepTolerance',1e-6,'ConstraintTolerance',1e-2, 'MaxIterations', 200);

sol_optimal=fmincon(objective,sol_0,[],[],[],[],lowbound,upbound,nonlcon,options);

theta_kick_optimal=sol_optimal(1)
h_pitchover_optimal=sol_optimal(2)
m_pld_max=sol_optimal(3)

function [c,ceq] = constraints(sol,stage,R,g,h_target,v_target,~,~)
theta_kick=sol(1);
h_pitchover=sol(2);
m_pld=sol(3);

[h_SECO, v_SECO, theta_SECO] = trajectory(theta_kick, h_pitchover, m_pld,stage,R,g);

    % equality constraints: hit target altitude, velocity, horizontal theta=0
    ceq = [ h_SECO - h_target;
            v_SECO - v_target;
            theta_SECO - 0; 
            ];
    c = []; % no inequality constraints
end 

%% simulator
function [h_SECO, v_SECO, theta_SECO] = trajectory(theta_kick, h_pitchover,m_pld,stage,R,g,v_orbit)

%% Burn 1

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

%h_pitchover=500; %height at start of gravity turn (m)
%theta_kick=deg2rad(1.3);

options=odeset('Events',@(t,Y) pitchover(t,Y,h_pitchover),'MaxStep',0.5);

Y0_1=[0;0;m0;0;pi/2]; %h=0, v=0, m=m0, x=0, theta=90
[t1,Y1,te,Ye,ie]=ode45(@(t,Y1) burn_1(t,Y1,stage,m_dot_1,R,g),[0 t_burn_1],Y0_1, options);

if ~isempty(ie)
    
    % Modify theta by applying the kick. To pitch over from vertical decrease theta:
    Ye(1,5) = Ye(1,5) - theta_kick;
    % Continue integration from event time to tfinal
    t_after = te(1);
    Y_after = Ye(1,:)';
    [t2, Y2] = ode45(@(t,Y) burn_1(t,Y,stage,m_dot_1,R,g), [t_after t_burn_1], ...
        Y_after, odeset(options,'InitialStep',1e-3));
    % concatenate results but avoid duplicating the event point
    t_1 = [t1(1:end-1); t2];
    Y_1 = [Y1(1:end-1,:); Y2];
else
    % no pitchover triggered (e.g., didn't reach h_pitchover during burn)
    t_1 = t1;
    Y_1 = Y1;
end


h_1=Y_1(:,1);
v_1=Y_1(:,2);
m_1=Y_1(:,3);
x_1=Y_1(:,4);
theta_1=Y_1(:,5);

v_MECO=v_1(end);
h_MECO=h_1(end);
m_MECO=stage(2).m_dry + stage(2).m_prop + m_pld;
x_MECO=x_1(end);
theta_MECO=theta_1(end);

function [value, isterminal, direction] = pitchover(~, Y, h_pitchover)
    % Event triggers when altitude crosses h_pitchover upward
    h = Y(1);
    value = h - h_pitchover; %event occurs when this equals 0
    isterminal = 1;    % stop integration 
    direction = +1;    % only when value increasing (crossing upward)
end

function dY1dt= burn_1(~,Y1,stage,m_dot_1,R,g)
h=Y1(1);
v=Y1(2);
m=Y1(3);
x=Y1(4);
theta=Y1(5);

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

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*(g*sin(theta)))/m;
dmdt=m_dot_1;
dthetadt=(-g/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY1dt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% stage seperation

t_seperation=7; %varies from 5-10 seconds
t_ignition_2=t_burn_1+t_seperation;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,R,g),[t_burn_1 t_ignition_2],Y0_sep);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
x_sep=Y_sep(:,4);
theta_sep=Y_sep(:,5);

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO;
x_SES=x_sep(end);
theta_SES=theta_sep(end);

function dY_sepdt= seperation(~,Y_sep,stage,R,g)
h=Y_sep(1);
v=Y_sep(2);
m=Y_sep(3);
x=Y_sep(4);
theta=Y_sep(5);

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

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g)/m;
dmdt=0;
dthetadt=(-g/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY_sepdt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% 2nd stage burn 

m_dot_2=-stage(2).thrust/(stage(2).isp * g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_SECO=t_burn_2+t_burn_1+t_seperation;

Y0_2=[h_SES;v_SES;m_SES;x_SES;theta_SES]; %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,R,g),[t_ignition_2 t_SECO],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);
x_2=Y2(:,4);
theta_2=Y2(:,5);

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
theta_SECO=theta_2(end);

function dY2dt= burn_2(~,Y2,stage,m_dot_2,R,g)
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

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g)/m;
dmdt=m_dot_2;
dthetadt=(-g/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY2dt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end


end 
