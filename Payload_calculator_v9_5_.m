
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

h_orbit=1200e3;
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

azimuth=deg2rad(90); %due east 

vehicle=Atlas_V_401();
stage=vehicle.stage;

burn_time_max=stage(2).m_prop/(stage(2).thrust/(stage(2).isp*g));
pld_int=stage(1).pld_guess+(stage(2).pld_guess-stage(1).pld_guess);

%% Optimisation

sol_0=   [ pld_int            ; 10 ; deg2rad(2)  ; 20 ; 7  ; 500          ; 2]; %initial guess for m_pld , t_seperation , theta_kick , t_pitchover ,  t_pitch
lowbound=[ stage(1).pld_guess ; 2  ; deg2rad(1)  ; 10 ; 5  ; 60           ; deg2rad(0)]; %best guess for lower bounds
upbound= [ stage(2).pld_guess ; 20 ; deg2rad(5)  ; 30 ; 10 ; burn_time_max; deg2rad(10)]; %best guess for upper bounds

objective=@(sol) -sol(1)/1e3;
nonlcon=@(sol) constraints(sol,stage,g,mu,h_orbit,r_earth,v_earth);

options=optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display','iter','StepTolerance',1e-8,'ConstraintTolerance',1e-4, 'MaxIterations', 1000);

problem=createOptimProblem('fmincon','x0',sol_0,'objective',objective,'lb',lowbound,'ub',upbound,'nonlcon',nonlcon,'options',options);
ms=MultiStart('UseParallel',false,'Display','iter');

num_starts=5;
rand_starts=rand(num_starts,length(sol_0)).*(upbound-lowbound)'+lowbound';
start_points=CustomStartPointSet(rand_starts);

[sol_optimal,fval,exitflag,output,solutions]=run(ms,problem,start_points);

m_pld_max=sol_optimal(1);
t_seperation_optimal=sol_optimal(2);
theta_kick_optimal=sol_optimal(3);
t_pitchover_optimal=sol_optimal(4);
t_pitch_optimal=sol_optimal(5);
t_burn_2a_optimal=sol_optimal(6);
theta_SECO_optimal=sol_optimal(7);

function [c,ceq] = constraints(sol,stage,g,mu,h_target,r_earth,v_earth)

m_pld=sol(1);
t_seperation=sol(2);
theta_kick=sol(3);
t_pitchover=sol(4);
t_pitch=sol(5);
t_burn_2a=sol(6);
theta_SECO=sol(7);

[h_SECO, v_SECO, m_SECO] = trajectory(m_pld,stage,g,mu,t_seperation,t_pitchover, theta_kick,t_pitch,t_burn_2a, theta_SECO);

%% Orbit calcs
r0 = [r_earth + h_SECO; 0];
v0 = [v_SECO*sin(theta_SECO); v_SECO*cos(theta_SECO) + v_earth]; % inertial

[h_SECO, v_SECO, m_SECO] = trajectory(m_pld,stage,g,mu,t_seperation,t_pitchover, theta_kick,t_pitch,t_burn_2a, theta_SECO);

%% Orbit calcs
r0 = [r_earth + h_SECO; 0];
v0 = [v_SECO*sin(theta_SECO); v_SECO*cos(theta_SECO) + v_earth]; % inertial

[dt_ap, r_ap, v_ap] = coast(r0, v0, mu);

h_ap=norm(r_ap)-r_earth;
r_ap_mag = norm(r_ap);
er = r_ap/r_ap_mag;
et = [-er(2); er(1)];

v_ap_r = dot(v_ap, er);   % should be ~0 at apogee
v_ap_t = abs(dot(v_ap, et));   % tangential speed at apogee

v_circ = sqrt(mu/r_ap_mag);
delta_v_2b = max(v_circ - v_ap_t, 0);
m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage(2).isp*g)));
m_prop_remaining=m_SECO-m_pld-stage(2).m_dry;

r0mag = norm(r0);
v0mag = norm(v0);

eps = v0mag^2/2 - mu/r0mag;          % specific energy
hvec = r0(1)*v0(2) - r0(2)*v0(1);    % scalar angular momentum (2D cross)
h_abs    = abs(hvec);

e = sqrt(1 + 2*eps*h_abs^2/mu^2);
a = -mu/(2*eps);

rp = a*(1 - e);     % perigee radius
hp = rp - r_earth;  % perigee altitude
h_min=120e3;

h_cut_min = 150e3;           % choose 150–200 km
c_hcut = (h_cut_min - h_SECO)/1000;  % <=0 means h_SECO >= h_cut_min

    % equality constraints: hit target altitude, mass
    ceq =  [(h_ap - h_target)/1000;
           ];
    c = [(m_prop_2b-m_prop_remaining)/10;
         (h_min-hp)/1000;
         (h_cut_min - h_SECO)/1000
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

function [h_SECO,v_SECO,m_SECO] = trajectory( m_pld,stage,g,mu,t_seperation,t_pitchover, theta_kick,t_pitch,t_burn_2a,theta_SECO)

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
m_MECO=stage(2).m_dry + stage(2).m_prop + m_pld;
x_MECO=x_1b(end);
theta_MECO=theta_1b(end);


%% Optimiser stage seperation

t_ignition_2a=t_burn_1+t_seperation;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[~,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g),[t_burn_1 t_ignition_2a],Y0_sep);

h_sep=Y_sep(:,1);
v_sep=Y_sep(:,2);
x_sep=Y_sep(:,4);
theta_sep=Y_sep(:,5);

v_SES=v_sep(end);
h_SES=h_sep(end);
m_SES=m_MECO;
x_SES=x_sep(end);
theta_SES=theta_sep(end);

%% Optimiser 2nd stage burn 

m_dot_2a= -stage(2).thrust/(stage(2).isp*g);

t_SECO=t_burn_2a+t_burn_1+t_seperation;

dthetadt_SES=(-(g*((6371e3/(6371e3+h_SES))^2))/(v_SES+0.1)-((v_SES+0.1)/(h_SES+6371e3)))*cos(theta_SES);

Y0_2a=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t_2a,Y2a]=ode45(@(t,Y2a) burn_2(t,Y2a,stage,m_dot_2a,g,t_ignition_2a,theta_SES,dthetadt_SES,theta_SECO),[t_ignition_2a t_SECO],Y0_2a);

h_2a=Y2a(:,1);
v_2a=Y2a(:,2);
m_2a=Y2a(:,3);
x_2a=Y2a(:,4);

t_rel_2a = t_2a - t_ignition_2a;
A_2a = tan(theta_SES);
B_2a = dthetadt_SES * (1 + A_2a^2);
k_2a=-A_2a/B_2a;
theta_2a = atan( A_2a + B_2a * k_2a *t_rel_2a );

v_SECO=v_2a(end);
h_SECO=h_2a(end);
x_SECO=x_2a(end);
m_SECO=m_2a(end);


%% Optimiser Trajectory ODEs

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

function dY2adt= burn_2(t,Y2a,stage,m_dot_2a,g,t_ignition_2a,theta_SES,dthetadt_SES,theta_SECO)
h=Y2a(1);
v=Y2a(2);
m=Y2a(3);
x=Y2a(4);

%t_rel=t-t_burn_1;
T=stage(2).thrust;
%m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);

t_rel = t - t_ignition_2a;
a = tan(theta_SES);
b = dthetadt_SES * (1 + a^2);
k=(tan(theta_SECO)-a)/b;
theta = atan( a + b * k*t_rel );
theta = min(theta, max(theta_SES,0));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_2a;

dY2adt = [dhdt;dvdt;dmdt;dxdt];

end


end

function [t_ap, r_ap, v_ap] = coast(r0, v0, mu)
% Find time to apogee by finding when vr changes sign.

% r0, v0 are 2D inertial vectors in the orbital plane.

r0mag = norm(r0);
er0 = r0/r0mag;
vr0 = dot(v0, er0);

% If already descending or nearly zero radial velocity, apogee is now.
if vr0 <= 1e-6
    t_ap = 0;
    r_ap = r0;
    v_ap = v0;
    return;
end

% Bracket the sign change (vr goes + to -)
t_lo = 0;
t_hi = 2000; % seconds (adjust if needed; LEO apogee is usually <~ 2000 s away)

for k = 1:30
    [r_hi, v_hi] = kepler_universal_2d(r0, v0, t_hi, mu);
    er_hi = r_hi/norm(r_hi);
    vr_hi = dot(v_hi, er_hi);
    if vr_hi < 0
        break;
    end
    t_hi = t_hi*1.5;
end

% Bisection to solve vr(t)=0
for k = 1:50
    t_mid = 0.5*(t_lo+t_hi);
    [r_mid, v_mid] = kepler_universal_2d(r0, v0, t_mid, mu);
    er_mid = r_mid/norm(r_mid);
    vr_mid = dot(v_mid, er_mid);

    if vr_mid > 0
        t_lo = t_mid;
    else
        t_hi = t_mid;
    end

    if abs(t_hi - t_lo) < 1e-3
        break;
    end
end

t_ap = 0.5*(t_lo+t_hi);
[r_ap, v_ap] = kepler_universal_2d(r0, v0, t_ap, mu);
end

function [r, v] = kepler_universal_2d(r0, v0, dt, mu)
% Universal variable Kepler propagation (works for elliptic/hyperbolic).

% 2D version using the same math as 3D.

r0mag = norm(r0);
v0sq  = dot(v0, v0);
vr0   = dot(r0, v0) / r0mag;

alpha = 2/r0mag - v0sq/mu;   % 1/a

% Initial guess for chi
if abs(alpha) > 1e-8
    chi = sqrt(mu)*abs(alpha)*dt;
else
    % near-parabolic
    chi = sqrt(mu)*dt/r0mag;
end

% Newton solve for chi
for k = 1:50
    z = alpha*chi^2;
    [C,S] = stumpff_CS(z);

    F = (r0mag*vr0/sqrt(mu))*chi^2*C + (1 - alpha*r0mag)*chi^3*S + r0mag*chi - sqrt(mu)*dt;
    dF = (r0mag*vr0/sqrt(mu))*chi*(1 - z*S) + (1 - alpha*r0mag)*chi^2*C + r0mag;

    dchi = -F/dF;
    chi = chi + dchi;

    if abs(dchi) < 1e-10
        break;
    end
end

z = alpha*chi^2;
[C,S] = stumpff_CS(z);

f = 1 - (chi^2/r0mag)*C;
g = dt - (1/sqrt(mu))*chi^3*S;

r = f*r0 + g*v0;
rmag = norm(r);

fdot = (sqrt(mu)/(rmag*r0mag))*(z*S - 1)*chi;
gdot = 1 - (chi^2/rmag)*C;

v = fdot*r0 + gdot*v0;
end

function [C,S] = stumpff_CS(z)
% Stumpff functions C(z), S(z)
if z > 1e-8
    sz = sqrt(z);
    C = (1 - cos(sz))/z;
    S = (sz - sin(sz))/(sz^3);
elseif z < -1e-8
    sz = sqrt(-z);
    C = (cosh(sz) - 1)/(-z);
    S = (sinh(sz) - sz)/(sz^3);
else
    % series expansions near 0
    C = 1/2;
    S = 1/6;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Burn 1a (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld_max;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1= -stage(1).m_prop / m_dot_1;

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

%% concatenate burn 1

t1=[t1a; t1b];
h_1=[h_1a; h_1b];
v_1=[v_1a; v_1b];
m_1=[m_1a; m_1b];
x_1=[x_1a; x_1b];
theta_1=[theta_1a; theta_1b];

%% stage seperation

t_ignition_2a=t_burn_1+t_seperation_optimal;

Y0_sep=[h_MECO;v_MECO;m_MECO;x_MECO;theta_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,g),[t_burn_1 t_ignition_2a],Y0_sep);

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
  

%% Burn 2a (up to apogee) 

m_dot_2a= -stage(2).thrust/(stage(2).isp*g);

t_SECO_a=t_burn_2a_optimal+t_burn_1+t_seperation_optimal;

dthetadt_SES=(-(g*((6371e3/(6371e3+h_SES))^2))/(v_SES+0.1)-((v_SES+0.1)/(h_SES+6371e3)))*cos(theta_SES);

Y0_2_a=[h_SES;v_SES;m_SES;x_SES]; %state at SES

[t_2a,Y2_a]=ode45(@(t,Y2_a) burn_2(t,Y2_a,stage,m_dot_2a,g,t_ignition_2a,theta_SES,dthetadt_SES,theta_SECO_optimal),[t_ignition_2a t_SECO_a],Y0_2_a);

h_2a=Y2_a(:,1);
v_2a=Y2_a(:,2);
m_2a=Y2_a(:,3);
x_2a=Y2_a(:,4);

t_rel_2a = t_2a - t_ignition_2a;
A_2a = tan(theta_SES);
B_2a = dthetadt_SES * (1 + A_2a^2);
k_2a=(tan(theta_SECO_optimal)-A_2a)/B_2a;
theta_2a = atan( A_2a + B_2a *k_2a*t_rel_2a );

v_SECO=v_2a(end);
h_SECO=h_2a(end);
x_SECO=x_2a(end);
m_SECO=m_2a(end);
theta_SECO=theta_2a(end);

%% redo orbit calcs
r0 = [r_earth + h_SECO; 0];
v0 = [v_SECO*sin(theta_SECO); v_SECO*cos(theta_SECO) + v_earth]; % inertial

[dt_ap, r_ap, v_ap] = coast(r0, v0, mu);

h_ap=norm(r_ap)-r_earth;
r_ap_mag = norm(r_ap);
er = r_ap/r_ap_mag;
et = [-er(2); er(1)];

v_ap_r = dot(v_ap, er);   % should be ~0 at apogee
v_ap_t = dot(v_ap, et);   % tangential speed at apogee

v_circ = sqrt(mu/r_ap_mag);
delta_v_2b = max(v_circ - v_ap_t, 0);
m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage(2).isp*g)));
m_prop_remaining=m_SECO-m_pld_max-stage(2).m_dry;

r0mag = norm(r0);
v0mag = norm(v0);
eps = v0mag^2/2 - mu/r0mag;          % specific energy
hvec = r0(1)*v0(2) - r0(2)*v0(1);    % scalar angular momentum (2D cross)
h_abs    = abs(hvec);
e = sqrt(1 + 2*eps*h_abs^2/mu^2);
a = -mu/(2*eps);
rp = a*(1 - e);     % perigee radius
hp = rp - r_earth;  % perigee altitude

%% ODE functions

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

T=stage(1).thrust-(P_atm*(stage(1).A_exit));

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

function dY2adt= burn_2(t,Y2a,stage,m_dot_2a,g,t_ignition_2a,theta_SES,dthetadt_SES,theta_SECO_optimal)
h=Y2a(1);
v=Y2a(2);
m=Y2a(3);
x=Y2a(4);

%t_rel=t-t_burn_1;
T=stage(2).thrust;
%m=stage(2).m_dry + stage(2).m_prop + m_pld + m_dot_2*t_rel;
Cd=stage(2).Cd;
A=stage(2).A;

[rho_atm,~]=isa_model(h);

t_rel = t - t_ignition_2a;
a = tan(theta_SES);
b = dthetadt_SES * (1 + a^2);
k=(tan(theta_SECO_optimal)-a)/b;
theta = atan( a + b * k*t_rel );
theta = min(theta, max(theta_SES,0));

dhdt=v*sin(theta);
dxdt=v*cos(theta);
dvdt=(T-(0.5*rho_atm*v^2*Cd*A)-m*g*((6371e3/(6371e3+h))^2)*sin(theta))/m;
dmdt=m_dot_2a;

dY2adt = [dhdt;dvdt;dmdt;dxdt];

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
theta_kick_optimal_deg=rad2deg(theta_kick_optimal);

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

fprintf('t_burn_2 = %.2f s\n', t_burn_2a_optimal)
fprintf('theta_SECO = %.2f degrees\n', rad2deg(theta_SECO))
fprintf('max_q = %.2f kPa\n', max_q)

fprintf('t_seperation_optimal = %.2f s\n', t_seperation_optimal)
fprintf('theta_kick_optimal = %.2f degs\n', theta_kick_optimal_deg)
fprintf('t_pitchover_optimal = %.2f s\n', t_pitchover_optimal)
fprintf('t_pitch_optimal = %.2f s\n', t_pitch_optimal)

fprintf('m_prop_2b = %.2f kg\n', m_prop_2b)
fprintf('m_prop_remaining = %.2f kg\n', m_prop_remaining)
fprintf('dt_ap = %.2f s\n', dt_ap)
fprintf('h_ap = %.2f km\n', h_ap/1000)
fprintf('h_per = %.2f km\n', hp/1000)
fprintf('e = %.2f \n', e)
fprintf('m_pld_max = %.2f kg\n', m_pld_max+(m_prop_remaining-m_prop_2b))