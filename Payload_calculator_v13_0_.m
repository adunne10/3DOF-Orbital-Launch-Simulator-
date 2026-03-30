
close all
clc
clear

%% Constants
global g G m_earth r_earth mu h_orbit v_earth vehicle t_seperation omega_earth

g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m
mu=G*m_earth;
omega_earth=7.2921159e-5; %units=rads/s

%% Mission Setup

vehicle=Atlas_V_401();
lat_launch=deg2rad(28); %cape canaveral
azimuth=deg2rad(90);

%% Loop
altitude=vehicle.performance.altitude;
interval=vehicle.performance.interval;
payload=vehicle.performance.payload;
h_vec = altitude(1):interval:altitude(end);     % meters
stage=vehicle.stage;

N = numel(h_vec);
m_pld_max_vec = nan(N,1);       % payload result (kg)
h_apo_vec     = nan(N,1);       % achieved apogee (m)
h_per_vec     = nan(N,1);       % achieved perigee (m)
e_vec         = nan(N,1);       % eccentricity
data_all(1:N) = struct('t',[],'x',[],'z',[],'vx',[],'vz',[],'m',[],'h',[],'theta',[],'v_inf',[]);
sols          =nan(N,6);
reserves      =nan(N,1);

for k = 1:N
    
%% Setup

%specify mission 
h_orbit = h_vec(k);

azimuth=deg2rad(90); %due east 

v_earth=omega_earth*r_earth*cos(lat_launch); %pure eastward velocity
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

t_seperation=10;

burn_time_max=stage(2).m_prop/(stage(2).thrust/(stage(2).isp*g));
burn_time_guess=300+((burn_time_max-300)/2);
pld_guess=stage(1).pld_guess+((stage(2).pld_guess-stage(1).pld_guess)/2);


%% Optimisation

%scale variables to around 1 
sol_scale = [pld_guess, 0.1, 10, 500, 10, 0.1];

%initial guesses
sol_0=   [ 
    pld_guess; %m_pld
    deg2rad(1); %theta_kick
    10; %t_pitchover
    burn_time_guess; %burn_time_2
    10; %t_pitch
    deg2rad(0) %alpha
]./sol_scale';

lowbound=[ 
    stage(1).pld_guess; %m_pld
    deg2rad(0); %theta_kick
    5; %t_pitchover
    200; %burn_time_2
    5; %t_pitch
    deg2rad(-6) %alpha
]./sol_scale';

upbound= [ 
    stage(2).pld_guess; %m_pld
    deg2rad(5); %theta_kick
    15;  %t_pitchover
    burn_time_max;  %burn_time_2
    15;  %t_pitch
    deg2rad(6); %alpha
]./sol_scale'; 

lambda = 0.8;  % trade off parameter
objective = @(sol) obj_payload_leftover(sol, sol_scale, h_orbit, lambda,lat_launch,azimuth);

nonlcon=@(sol) constraints(sol,sol_scale,h_orbit,lat_launch,azimuth);

options=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'FiniteDifferenceType','forward' , ...
    'Display','iter', ...
    'FiniteDifferenceStepSize',1e-3, ...
    'ConstraintTolerance',1e-2, ...
    'OptimalityTolerance', 1e-4, ...
    'MaxIterations', 2000 ...
    );

% if k > 1 && all(isfinite(sols(k-1,:)))
%     sol_0 = (sols(k-1,:)' ./ sol_scale');   % previous solution, scaled
% end

sol__optimal_scaled=fmincon(objective,sol_0,[],[],[],[],lowbound,upbound,nonlcon,options);

sol_optimal=sol__optimal_scaled.*sol_scale';

%% Results 

data=sim(sol_optimal,lat_launch,azimuth);
[h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit,sol_optimal);

%plot_orbit_from_SECO(sol_optimal, data, h_orbit);

m_pld_max = sol_optimal(1);
m_reserve=(m_prop_remaining - m_prop_2b);

% store
m_pld_max_vec(k) = m_pld_max;
h_apo_vec(k)     = h_apo;
h_per_vec(k)     = h_per;
e_vec(k)         = e;
data_all(k)      =data;
sols(k,:)        =sol_optimal(:).';
reserves(k)      =m_reserve;

end

error=100.*((m_pld_max_vec-payload)./payload);

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
    plot(data_all(k).t,(data_all(k).h./1000),'DisplayName', sprintf('h=%.0f km', h_vec(k)/1000));
end

% yline(pi()/2)
xlabel('Time (s)')
ylabel('altitude (km)')
legend('Location','best')
hold off

figure; hold on; grid on

for k = 1:N
    plot(data_all(k).t,(data_all(k).v_inf),'DisplayName', sprintf('h=%.0f km', h_vec(k)/1000));
end

% yline(pi()/2)
xlabel('Time (s)')
ylabel('speed (m/s)')
legend('Location','best')
hold off


figure;
bar(h_vec'./1000,error)
grid
xlabel('Target orbit altitude (km)')
ylabel('error (%)')
hold on
yline(10,'r');
yline(-10,'r');
hold off

%% Script functions

function [c,ceq] = constraints(sol,sol_scale,h_orbit,lat_launch,azimuth)

sol_unscaled=sol.*sol_scale';

data =sim(sol_unscaled,lat_launch,azimuth);

[h_apo, h_per, ~, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit, sol_unscaled);

h_min = 120e3;
h_cut_min = 150e3;         
h_SECO=data.h(end);
theta_SECO=data.theta(end);
theta_SECO_max=deg2rad(5);
theta_SECO_min=deg2rad(-5);


    % equality constraints: hit target altitude, mass
    ceq =  [(h_apo - h_orbit)/1000;
           ];

    c = [(m_prop_2b-m_prop_remaining);
         (h_min-h_per)/1000;
         (h_cut_min - h_SECO)/1000;
         theta_SECO_min - theta_SECO;   
         theta_SECO - theta_SECO_max    
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

function data =sim(sol,lat_launch,azimuth)

global g vehicle t_seperation r_earth omega_earth

stage=vehicle.stage;

m_pld=sol(1);
theta_kick=sol(2);
t_pitchover=sol(3);
t_burn_2=sol(4);
t_pitch=sol(5);
alpha=sol(6);

opts_ode = odeset('RelTol',1e-8,'AbsTol',1e-8);
   
% Burn 1a (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1= -stage(1).m_prop / m_dot_1;
t_burn_1a=t_pitchover+t_pitch;
vx0 = omega_earth * r_earth * cos(lat_launch) * sin(azimuth);

Y0_1a=[0; r_earth; vx0; 0; m0]; %x=0, z=r_earth, vx=depends on location, vz=0, m=GLOW

[t1a,Y1a]=ode45(@(t,Y1a) burn_1a(t,Y1a,stage,m_dot_1,t_pitchover,theta_kick,t_pitch,lat_launch,azimuth),[0 t_burn_1a],Y0_1a,opts_ode);

x_1a=Y1a(:,1); 
z_1a=Y1a(:,2);
vx_1a=Y1a(:,3); 
vz_1a=Y1a(:,4);

r_1a=hypot(x_1a,z_1a);
urx_1a=x_1a./r_1a;  
urz_1a=z_1a./r_1a;
utx_1a=urz_1a;   
utz_1a=-urx_1a;

theta_1a = nan(size(t1a));

for i = 1:numel(t1a)
    if t1a(i) < t_pitchover
        % thrust straight up = ur
        tx_1a = urx_1a(i);  tz_1a = urz_1a(i);
    elseif t1a(i) < t_pitchover + t_pitch
        ramp  = (t1a(i) - t_pitchover)/t_pitch;
        theta = ramp * theta_kick;
        tx_1a = cos(theta)*urx_1a(i) + sin(theta)*utx_1a(i);
        tz_1a = cos(theta)*urz_1a(i) + sin(theta)*utz_1a(i);
    else
        % thrust aligned with velocity
        vmag_1a = hypot(vx_1a(i), vz_1a(i));
        tx_1a = vx_1a(i)/vmag_1a;  tz_1a = vz_1a(i)/vmag_1a;
    end

    % project thrust onto local (ur, ut)
    tr_1a = tx_1a*urx_1a(i) + tz_1a*urz_1a(i);
    tt_1a = tx_1a*utx_1a(i) + tz_1a*utz_1a(i);

    theta_1a(i) = atan2(tr_1a, tt_1a);   % 0=horizontal, pi/2=up
end

% Burn 1b (gravity turn)

t_burn_1b=t_burn_1-t_burn_1a;
Y0_1b=Y1a(end,:)'; 

[t1b,Y1b]=ode45(@(t,Y1b) burn_1b(t,Y1b,stage,m_dot_1,lat_launch,azimuth),[t_burn_1a t_burn_1b+t_burn_1a],Y0_1b,opts_ode);

x_1b=Y1b(:,1);
z_1b=Y1b(:,2);
vx_1b=Y1b(:,3);
vz_1b=Y1b(:,4);
m_1b=Y1b(:,5);

r_1b= hypot(x_1b,z_1b);
urx_1b=x_1b./r_1b;  
urz_1b=z_1b./r_1b;
utx_1b=urz_1b;   
utz_1b=-urx_1b;

vmag_1b = hypot(vx_1b, vz_1b);
tx_1b = vx_1b./vmag_1b;
tz_1b = vz_1b./vmag_1b;

tr_1b = tx_1b.*urx_1b + tz_1b.*urz_1b;
tt_1b = tx_1b.*utx_1b + tz_1b.*utz_1b;

theta_1b = atan2(tr_1b, tt_1b);

x_MECO = x_1b(end);
z_MECO = z_1b(end);
vx_MECO = vx_1b(end);
vz_MECO = vz_1b(end);
m_MECO = m_1b(end) - stage(1).m_dry;

% stage seperation

t_ignition_2=t_burn_1+t_seperation;

Y0_sep=[x_MECO;z_MECO;vx_MECO;vz_MECO;m_MECO]; %state at MECO 

[t_sep,Y_sep]=ode45(@(t,Y_sep) seperation(t,Y_sep,stage,lat_launch,azimuth),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

x_sep = Y_sep(:,1); 
z_sep = Y_sep(:,2);
vx_sep = Y_sep(:,3); 
vz_sep = Y_sep(:,4);

r_sep  = hypot(x_sep,z_sep);
urx_sep = x_sep./r_sep;  
urz_sep = z_sep./r_sep;
utx_sep = urz_sep;   
utz_sep = -urx_sep;

vmag_sep = hypot(vx_sep,vz_sep);
tx_sep = vx_sep./vmag_sep;  
tz_sep = vz_sep./vmag_sep;

tr_sep = tx_sep.*urx_sep + tz_sep.*urz_sep;
tt_sep = tx_sep.*utx_sep + tz_sep.*utz_sep;

theta_sep = atan2(tr_sep, tt_sep);

% Burn 2 (up to apogee) 

m_dot_2= -stage(2).thrust/(stage(2).isp*g);
t_SECO=t_burn_2+t_burn_1+t_seperation;

Y0_2=Y_sep(end,:)';  %state at SES

[t2,Y2]=ode45(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,lat_launch,azimuth,alpha),[t_ignition_2 t_SECO],Y0_2,opts_ode);

x_2 = Y2(:,1); 
z_2 = Y2(:,2);
vx_2 = Y2(:,3); 
vz_2 = Y2(:,4);

r_2  = hypot(x_2,z_2);
urx_2 = x_2./r_2;  
urz_2 = z_2./r_2;
utx_2 = urz_2;   
utz_2 = -urx_2;

vmag_2 = hypot(vx_2,vz_2);
ux_2 = vx_2./vmag_2;  
uz_2 = vz_2./vmag_2;
nx_2 = -uz_2;        
nz_2 = ux_2;

tx_2 = cos(alpha).*ux_2 + sin(alpha).*nx_2;
tz_2 = cos(alpha).*uz_2 + sin(alpha).*nz_2;

tr_2 = tx_2.*urx_2 + tz_2.*urz_2;
tt_2 = tx_2.*utx_2 + tz_2.*utz_2;

theta_2 = atan2(tr_2, tt_2);

% Data output

data.t=     [t1a;       t1b(2:end);      t_sep(2:end);      t2(2:end)];
data.x=     [Y1a(:,1);  Y1b(2:end,1);    Y_sep(2:end,1);    Y2(2:end,1)];
data.z=     [Y1a(:,2);  Y1b(2:end,2);    Y_sep(2:end,2);    Y2(2:end,2)];
data.vx=    [Y1a(:,3);  Y1b(2:end,3);    Y_sep(2:end,3);    Y2(2:end,3)];
data.vz=    [Y1a(:,4);  Y1b(2:end,4);    Y_sep(2:end,4);    Y2(2:end,4)];
data.m =    [Y1a(:,5);  Y1b(2:end,5);    Y_sep(2:end,5);    Y2(2:end,5)];
data.theta= [theta_1a;  theta_1b(2:end); theta_sep(2:end);  theta_2(2:end)];      % flight-path angle (0 = horizontal)
data.h=     sqrt(data.x.^2 + data.z.^2)-r_earth;     
    r_1b = sqrt(data.x.^2 + data.z.^2); 
    v_atm = omega_earth .* r_1b .* cos(lat_launch) .* sin(azimuth);

data.v_inf= sqrt((data.vx-v_atm).^2 + data.vz.^2);

end 

function [h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, ~ ,sol)

 global mu r_earth g vehicle

stage=vehicle.stage;
m_pld=sol(1);

x_SECO=data.x(end);
z_SECO=data.z(end);
vx_SECO=data.vx(end);
vz_SECO=data.vz(end);
m_SECO=data.m(end);

r0 = [x_SECO; z_SECO];
v0 = [vx_SECO; vz_SECO]; % inertial

r0mag = norm(r0);
v0mag = norm(v0);

orbit_energy_sp  = v0mag^2/2 - mu/r0mag;
h_vec = r0(1)*v0(2) - r0(2)*v0(1);
h_abs = abs(h_vec);

e = sqrt(1 + 2*orbit_energy_sp*h_abs^2/mu^2);
a = -mu/(2*orbit_energy_sp);

r_per = a*(1 - e);
r_apo = a*(1 + e);

h_per = r_per - r_earth;
h_apo = r_apo - r_earth;

% circularization at apogee 
va = sqrt(mu*(2/r_apo - 1/a));
vcirc = sqrt(mu/r_apo);
delta_v_2b = max(vcirc - va, 0);

m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage(2).isp*g)));
m_prop_remaining=m_SECO-m_pld-stage(2).m_dry;

end

function J = obj_payload_leftover(sol, sol_scale, h_orbit, lambda,lat_launch,azimuth)

sol_unscaled = sol .* sol_scale';

data = sim(sol_unscaled,lat_launch,azimuth);
[~, ~, ~, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit, sol_unscaled);

m_pld  = sol_unscaled(1);
m_left = max(m_prop_remaining - m_prop_2b, 0);  % leftover after planned circ (kg)

% maximize (m_pld - lambda*m_left)  <=>  minimize negative:
J = -(m_pld - lambda*m_left);

end

function plot_orbit_from_SECO(~, data, h_orbit)
    % Visualize the 2-body orbit implied by SECO state, plus target circle.
    global mu r_earth v_earth

    h_SECO     = data.h(end);
    v_SECO     = data.v(end);
    theta_SECO = data.theta(end);

    % Inertial position/velocity at SECO (same mapping you use in orbit_calc)
    r0 = [r_earth + h_SECO; 0];
    v0 = [v_SECO*sin(theta_SECO); v_SECO*cos(theta_SECO) + v_earth];

    r0mag = norm(r0);
    v0mag = norm(v0);

    % Specific angular momentum (scalar in 2D)
    h = r0(1)*v0(2) - r0(2)*v0(1);

    % Specific orbital energy
    eps = v0mag^2/2 - mu/r0mag;

    % Orbit elements
    e = sqrt(1 + 2*eps*h^2/mu^2);
    a = -mu/(2*eps);
    p = h^2/mu;

    % True anomaly at r0: use polar form r = p/(1 + e cos nu)
    cosnu0 = (p/r0mag - 1)/e;
    cosnu0 = max(min(cosnu0,1),-1); % clamp
    nu0 = acos(cosnu0);
    % Choose sign based on radial velocity (outbound/inbound)
    vr0 = dot(r0, v0)/r0mag;
    if vr0 < 0, nu0 = -nu0; end

    % Eccentricity vector to get periapsis direction
    evec = ( (v0mag^2 - mu/r0mag).*r0 - dot(r0,v0).*v0 )/mu;
    omega = atan2(evec(2), evec(1)); % argument of periapsis in inertial xy

    % Sample true anomaly around the orbit
    nu = linspace(-pi, pi, 2000);
    r  = p ./ (1 + e*cos(nu));

    % Rotate perifocal->inertial by omega
    x = r .* cos(nu + omega);
    y = r .* sin(nu + omega);

    % Earth circle
    ang = linspace(0,2*pi,600);
    xe = r_earth*cos(ang);
    ye = r_earth*sin(ang);

    % Target circular orbit radius
    rt = r_earth + h_orbit;
    xt = rt*cos(ang);
    yt = rt*sin(ang);

    % Apogee / perigee points
    r_per = a*(1 - e);
    r_apo = a*(1 + e);
    xper = r_per*cos(omega);
    yper = r_per*sin(omega);
    xapo = -r_apo*cos(omega); % apo is opposite peri direction
    yapo = -r_apo*sin(omega);

    figure; hold on; axis equal; grid on
    plot(xe/1000, ye/1000, 'b', 'LineWidth', 1.5);          % Earth
    plot(xt/1000, yt/1000, '--', 'LineWidth', 1.2);         % target circle
    plot(x/1000,  y/1000,  'LineWidth', 1.5);               % transfer ellipse
    plot(r0(1)/1000, r0(2)/1000, '*', 'MarkerSize', 7);     % SECO point
    plot(xper/1000, yper/1000, '*', 'MarkerSize', 7);       % perigee
    plot(xapo/1000, yapo/1000, '*', 'MarkerSize', 7);       % apogee

    xlabel('Inertial x (km)')
    ylabel('Inertial y (km)')
    title(sprintf('Orbit from SECO: e=%.4f, a=%.0f km, h_tgt=%.0f km', e, a/1000, h_orbit/1000))
    legend('Earth','Target circular orbit','Transfer orbit','SECO','Perigee','Apogee','Location','best')
end


%% ODE functions

function dY1adt= burn_1a(t,Y1a,stage,m_dot_1,t_pitchover,theta_kick,t_pitch,lat_launch,azimuth)
x=Y1a(1);
z=Y1a(2);
vx=Y1a(3);
vz=Y1a(4);
m=Y1a(5);

global mu r_earth omega_earth 

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_inf=max(sqrt((vx-v_atm)^2 + vz^2),0.1);

[rho_atm,P_atm]=isa_model(h);
Cd=stage(1).Cd;
A=stage(1).A;

g_x=-(mu*x/r^3);
g_z=-(mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

urx = x/r;  %unit vector in radial 'up' direction (local vertical)
urz = z/r; 
utx = urz; % unit vector in tangental direction (local horizontal)
utz = -urx; 

if t < t_pitchover
    tx = urx; 
    tz = urz;           % straight up (radial)
elseif t < t_pitchover + t_pitch
    ramp = (t - t_pitchover)/t_pitch; 
    theta = ramp * theta_kick;          % ramp to theta_kick
    tx = cos(theta)*urx + sin(theta)*utx;
    tz = cos(theta)*urz + sin(theta)*utz;
else
    % after pitchover, align thrust with velocity 
    tx =vx/sqrt(vx^2 + vz^2); 
    tz =vz/sqrt(vx^2 + vz^2);
end

T=stage(1).thrust-(P_atm*(stage(1).A_exit));
T_x=T*tx;
T_z=T*tz;


dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_1;

dY1adt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY1bdt= burn_1b(~,Y1b,stage,m_dot_1,lat_launch,azimuth)

x=Y1b(1);
z=Y1b(2);
vx=Y1b(3);
vz=Y1b(4);
m=Y1b(5);

global mu r_earth omega_earth 

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_inf=max(sqrt((vx-v_atm)^2 + vz^2),0.1);

[rho_atm,P_atm]=isa_model(h);
Cd=stage(1).Cd;
A=stage(1).A;

g_x=-(mu*x/r^3);
g_z=-(mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

tx =vx/sqrt(vx^2 + vz^2); 
tz =vz/sqrt(vx^2 + vz^2);
T=stage(1).thrust-(P_atm*(stage(1).A_exit));
T_x=T*tx;
T_z=T*tz;

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_1;

dY1bdt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY_sepdt= seperation(~,Y_sep,stage,lat_launch,azimuth)
x=Y_sep(1);
z=Y_sep(2);
vx=Y_sep(3);
vz=Y_sep(4);
m=Y_sep(5);

global mu r_earth omega_earth 

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_inf=max(sqrt((vx-v_atm)^2 + vz^2),0.1);

[rho_atm,~]=isa_model(h);
Cd=stage(2).Cd;
A=stage(2).A;

g_x=-(mu*x/r^3);
g_z=-(mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(-D_x)/m;
dvzdt=g_z+(-D_z)/m;
dmdt=0;

dY_sepdt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY2dt= burn_2(~,Y2,stage,m_dot_2,lat_launch,azimuth,alpha)
x=Y2(1);
z=Y2(2);
vx=Y2(3);
vz=Y2(4);
m=Y2(5);

global mu r_earth omega_earth 

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_inf=max(sqrt((vx-v_atm)^2 + vz^2),0.1);

[rho_atm,P_atm]=isa_model(h);
Cd=stage(2).Cd;
A=stage(2).A;

g_x=-(mu*x/r^3);
g_z=-(mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

ux = vx/sqrt(vx^2 + vz^2); uz = vz/sqrt(vx^2 + vz^2);
nx = -uz;  nz = ux;

tx=cos(alpha)*ux + sin(alpha)*nx;
tz=cos(alpha)*uz + sin(alpha)*nz;
T=stage(2).thrust-(P_atm*(stage(2).A_exit));
T_x=T*tx;
T_z=T*tz;

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_2;

dY2dt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

    end
