
close all
clc
clear

%% Setup

R=287; %units=J/kg*K
g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m

h_orbit=650e3;
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

theta_SECO_goal=deg2rad(0.5);

stage=[
    struct('m_dry', 25000, 'm_prop', 400000, 'isp', 283, 'thrust', 7.6e6, ...
    'Cd', 0.3, 'A', 10, 'A_exit', 0.9)
    struct('m_dry', 3900, 'm_prop', 92000, 'isp', 348, 'thrust', 981e3, ...
    'Cd', 0.25, 'A', 10, 'A_exit',0) 
    ];


%% Optimisation

sol_0=   [10000 ; 10  ; deg2rad(1) ; 500]; %initial guess for m_pld ,  t_coast & theta_kick

lowbound=[ 5000 ; 5 ; deg2rad(0.1) ; 100]; %best guess for lower bounds
upbound= [ 22000 ; 30 ; deg2rad(2) ; 1000]; %best guess for upper bounds

objective=@(sol) -sol(1)/1e3;

nonlcon=@(sol) constraints(sol,stage,R,g,h_orbit,v_orbit,m_earth,r_earth,theta_SECO_goal);
options=optimoptions('fmincon', 'Algorithm', 'interior-point', 'EnableFeasibilityMode',true, 'Display','iter','StepTolerance',1e-8,'ConstraintTolerance',1e-4, 'OptimalityTolerance', 1e-4, 'MaxIterations', 2000);

sol_optimal=fmincon(objective,sol_0,[],[],[],[],lowbound,upbound,nonlcon,options);

m_pld_max=sol_optimal(1);
t_coast_optimal=sol_optimal(2);
theta_kick_optimal=sol_optimal(3);
h_pitchover_optimal=sol_optimal(4);

function [c,ceq] = constraints(sol,stage,R,g,h_target,v_target,~,~,theta_SECO_goal)

m_pld=sol(1);
t_coast=sol(2);
theta_kick=sol(3);
h_pitchover=sol(4);

[h_SECO, v_SECO] = trajectory(m_pld,stage,R,g,t_coast,h_pitchover, theta_kick,theta_SECO_goal);

    % equality constraints: hit target altitude, velocity, horizontal theta=0
    ceq =  [
            (h_SECO - h_target)/1000;
            (v_SECO - v_target)/100 
           ];
    c = [];
end 

function [h_SECO, v_SECO] = trajectory( m_pld,stage,R,g,t_coast,h_pitchover, theta_kick,theta_SECO_goal)

%% Optimiser Burn 1

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

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
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY1dt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% Optimiser stage seperation

t_ignition_2=t_burn_1+t_coast;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[~,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,R,g),[t_burn_1 t_ignition_2],Y0_sep);

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
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2))/m;
dmdt=0;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);


dY_sepdt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% Optimiser 2nd stage burn 

m_dot_2=-stage(2).thrust/(stage(2).isp*g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_SECO=t_burn_2+t_burn_1+t_coast;

dthetadt_SES=(-(g*((6371e3/(6371e3+h_SES))^2))/(v_SES+0.1)-((v_SES+0.1)/(h_SES+6371e3)))*cos(theta_SES);

Y0_2=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,R,g,dthetadt_SES,theta_SES,t_ignition_2,t_burn_2,theta_SECO_goal),[t_ignition_2 t_SECO],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);
x_2=Y2(:,4);

tau_2 = (t2 - t_ignition_2) / t_burn_2;
tau_2 = min(max(tau_2,0),1);

a = tan(theta_SES);
k = tan(theta_SECO_goal) - a;

theta_2 = atan(a + k * tau_2);

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
theta_SECO=theta_2(end);

    function dY2dt= burn_2(t,Y2,stage,m_dot_2,R,g,dthetadt_SES,theta_SES,t_ignition_2 ,t_burn_2,theta_SECO_goal)
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

tau = (t - t_ignition_2) / t_burn_2;
tau = min(max(tau,0),1);
a = tan(theta_SES);
b = dthetadt_SES * (1 + a^2);
k = tan(theta_SECO_goal) - a;     % since tau ends at 1
theta = atan( a + k * tau );

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2))/m;
dmdt=m_dot_2;

dY2dt = [dhdt;dvdt;dmdt;dxdt];

end


end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Burn 1

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld_max;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1=-stage(1).m_prop / m_dot_1;

options=odeset('Events',@(t,Y) pitchover(t,Y,h_pitchover_optimal),'MaxStep',0.5);

Y0_1=[0;0;m0;0;pi/2]; %h=0, v=0, m=m0, x=0, theta=90
[t1,Y1,te,Ye,ie]=ode45(@(t,Y1) burn_1(t,Y1,stage,m_dot_1,R,g),[0 t_burn_1],Y0_1, options);

if ~isempty(ie)
    
    % Modify theta by applying the kick. To pitch over from vertical decrease theta:
    Ye(1,5) = Ye(1,5) - theta_kick_optimal;
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
m_MECO=stage(2).m_dry + stage(2).m_prop + m_pld_max;
x_MECO=x_1(end);
theta_MECO=theta_1(end);

function [value, isterminal, direction] = pitchover(~, Y, h_pitchover_optimal)
    % Event triggers when altitude crosses h_pitchover upward
    h = Y(1);
    value = h - h_pitchover_optimal; %event occurs when this equals 0
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
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_1;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);

dY1dt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% stage seperation

t_ignition_2=t_burn_1+t_coast_optimal;

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
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2))/m;
dmdt=0;
dthetadt=(-(g*((6371e3/(6371e3+h))^2))/(v+0.1)-((v+0.1)/(h+6371e3)))*cos(theta);


dY_sepdt = [dhdt;dvdt;dmdt;dxdt;dthetadt];

end

%% 2nd stage burn 

m_dot_2=-stage(2).thrust/(stage(2).isp*g);
t_burn_2=-stage(2).m_prop / m_dot_2;
t_SECO=t_burn_2+t_burn_1+t_coast_optimal;

dthetadt_SES=(-(g*((6371e3/(6371e3+h_SES))^2))/(v_SES+0.1)-((v_SES+0.1)/(h_SES+6371e3)))*cos(theta_SES);

Y0_2=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,R,g,dthetadt_SES,theta_SES,t_ignition_2,t_burn_2,theta_SECO_goal),[t_ignition_2 t_SECO],Y0_2);

h_2=Y2(:,1);
v_2=Y2(:,2);
m_2=Y2(:,3);
x_2=Y2(:,4);

tau_2 = (t2 - t_ignition_2) / t_burn_2;
tau_2 = min(max(tau_2,0),1);
a = tan(theta_SES);
k = tan(theta_SECO_goal) - a;
theta_2 = atan(a + k * tau_2);

v_SECO=v_2(end);
h_SECO=h_2(end);
x_SECO=x_2(end);
theta_SECO=theta_2(end);

function dY2dt= burn_2(t,Y2,stage,m_dot_2,R,g,dthetadt_SES,theta_SES,t_ignition_2,t_burn_2,theta_SECO_goal)
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

tau = (t - t_ignition_2) / t_burn_2;
tau = min(max(tau,0),1);
a = tan(theta_SES);
b = dthetadt_SES * (1 + a^2);
k = tan(theta_SECO_goal) - a;     % since tau ends at 1
theta = atan( a + k * tau );

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

plot(t_1,h_1_km) %altitude vs time
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
plot(t_1,v_1) % veolcity vs time
grid
hold on
plot (t_sep,v_sep)
plot (t2,v_2)
ylabel("velocity (m/s)")
xlabel("time (s)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off 


figure%4
plot(t_1,theta_1_deg) % theta vs time
grid
hold on
plot (t_sep,theta_sep_deg)
plot (t2,theta_2_deg)
xlabel("time (s)")
ylabel("theta (deg)")
legend("Burn 1","Stage seperation","Burn 2",Location="southeast")
hold off 

% figure %5
% plot(q_1,h_1_km) %altitude vs dynamic pressure
% grid
% hold on
% xlabel("dynamic pressure (kPa)")
% ylabel("altitude (km)")
% legend("Burn 1",Location="southeast")
% xlim([0 60]);
% ylim([0 30]);
% hold off 
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

fprintf('t_coast_optimal = %.2f s\n', t_coast_optimal)
fprintf('theta_kick_optimal = %.2f rads\n', theta_kick_optimal)
fprintf('h_pitchover_optimal = %.2f m\n', h_pitchover_optimal)
fprintf('m_pld_max = %.2f kg\n', m_pld_max)

