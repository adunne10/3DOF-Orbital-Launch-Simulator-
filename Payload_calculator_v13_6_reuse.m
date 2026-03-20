
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

%% LV setup
vehicle=Falcon9();
altitude=vehicle.performance.altitude;
h_vec = vehicle.performance.altitude;     % meters
stage=vehicle.stage;

N = numel(h_vec);

inc_deg = [40 60 80 105];          % inclinations to sweep (deg)
Ni    = numel(inc_deg);

m_pld_max_vec_i  = nan(N, Ni);      % simulated max payload (kg) vs altitude, per inclination
m_pld_max_vec = nan(N,1);       % payload result (kg)
% h_apo_vec     = nan(N,1);       % achieved apogee (m)
% h_per_vec     = nan(N,1);       % achieved perigee (m)
% e_vec         = nan(N,1);       % eccentricity
% data_all(1:N) = struct('t',[],'x',[],'z',[],'vx',[],'vz',[],'m',[],'h',[],'theta',[],'v_inf',[]);
% sols          =nan(N,7);
% reserves      =nan(N,1);
% exitflag_vec = nan(N,1);
pld_guess_vec = nan(N,1);

%% Mission Setup

lat_launch=deg2rad(28); 


i_target=deg2rad(inc_deg(28)); %must be greater than abs(lat_launch)

if i_target<(pi/2)

azimuth=asin(cos(i_target)/cos(lat_launch));

else 
    azimuth = pi - asin(cos(i_target)/cos(lat_launch));
end

%% Loop
for k = 1:N
    
%% Setup
 
h_orbit = h_vec(k);
v_earth=omega_earth*r_earth*cos(lat_launch); %pure eastward velocity
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

dV_orbit=sqrt(mu/(r_earth+h_orbit)); %m/s 
dV_grav_loss=0.6*sqrt((2*g*h_orbit)/(1+(h_orbit/r_earth))); %m/s (approx.)
dV_drag=100; %m/s (approx.)
v_earth=omega_earth*r_earth*cos(lat_launch); %pure eastward velocity
dV_req=dV_orbit+dV_drag+dV_grav_loss-v_earth;

m_wet=stage(1).m_dry+stage(1).m_prop+stage(2).m_dry+stage(2).m_prop;

m_wet_1=stage(1).m_dry+stage(1).m_prop;
m_wet_2=stage(2).m_dry+stage(2).m_prop;

dv_residual = @(m_pld) ...
    stage(1).isp * g * log((m_wet_1 + m_wet_2 + m_pld) / ...
                           (stage(1).m_dry + m_wet_2 + m_pld)) + ...
    stage(2).isp * g * log((m_wet_2 + m_pld) / ...
                           (stage(2).m_dry + m_pld)) ...
    - dV_req ;

pld_guess = fzero(dv_residual, 0.02 * m_wet);

t_seperation=10;

burn_time_max=stage(2).m_prop/(stage(2).thrust/(stage(2).isp*g));
burn_time_guess=300+((burn_time_max-300)/2);

%% Optimisation

%scale variables to around 1 
sol_scale = [pld_guess, deg2rad(5), 10, burn_time_guess, 10, deg2rad(6), deg2rad(0.05)];

%initial guesses
sol_0=   [ 
    0.66*pld_guess; %m_pld
    deg2rad(1); %theta_kick
    10; %t_pitchover
    burn_time_guess; %burn_time_2
    10; %t_pitch
    deg2rad(0); %alpha_0
    deg2rad(0) %alpha rate
]./sol_scale';

lowbound=[ 
    (0.1*pld_guess); %m_pld
    deg2rad(0.5); %theta_kick
    5; %t_pitchover
    250; %burn_time_2
    5; %t_pitch
    deg2rad(-6); %alpha_0
    deg2rad(-0.05) %alpha rate
]./sol_scale';

upbound= [ 
    (7*pld_guess); %m_pld
    deg2rad(5); %theta_kick
    30;  %t_pitchover
    burn_time_max;  %burn_time_2
    20;  %t_pitch
    deg2rad(6); %alpha_0
    deg2rad(0.05) %alpha rate
]./sol_scale'; 

lambda = 0.2;  % trade off parameter
objective = @(sol) obj_payload_leftover(sol, sol_scale, h_orbit, lambda,lat_launch,azimuth);

nonlcon=@(sol) constraints(sol,sol_scale,h_orbit,lat_launch,azimuth);

options=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'FiniteDifferenceType','central' , ...
    'Display','iter', ...
    'FiniteDifferenceStepSize',1e-4, ...
    'ConstraintTolerance',5e-2, ...
    'OptimalityTolerance', 1e-3, ...
    'MaxIterations', 500, ...
    'MaxFunctionEvaluations', 5e5 ...
    );

[sol__optimal_scaled, fval, exitflag, output]=fmincon(objective,sol_0,[],[],[],[],lowbound,upbound,nonlcon,options);

sol_optimal=sol__optimal_scaled.*sol_scale';

%% Results 

data=sim(sol_optimal,lat_launch,azimuth);
[h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit,sol_optimal);

m_pld_max = sol_optimal(1);
m_reserve=(m_prop_remaining - m_prop_2b);

% store
m_pld_max_vec(k) = m_pld_max;

end


%% Plots

h_vec_km=h_vec'./1000;

figure; hold on; grid on;
hSim = gobjects(Ni,1);
% plot simulated curves (one per inclination)
for j = 1:Ni
    hSim(j) = plot(h_vec_km, m_pld_max_vec_i(:,j), '-*', ...
        'DisplayName', sprintf('Sim %d°', inc_deg(j)));
end
  plot(h_vec_km, vehicle.performance.payload_40, '--', ...
    'Color', hSim(inc_deg==40).Color, 'DisplayName','Actual 40°');

plot(h_vec_km, vehicle.performance.payload_60, '--', ...
    'Color', hSim(inc_deg==60).Color, 'DisplayName','Actual 60°');

plot(h_vec_km, vehicle.performance.payload_80, '--', ...
    'Color', hSim(inc_deg==80).Color, 'DisplayName','Actual 80°');

plot(h_vec_km, vehicle.performance.payload_SSO, '--', ...
    'Color', hSim(inc_deg==105).Color, 'DisplayName','Actual 105° (SSO)');
xlabel('Target orbit altitude (km)');
ylabel('Max payload (kg)');
legend('Location','best');
hold off;

% figure;hold on;grid on;
% plot(h_vec_km(ok),error(ok),'-x','MarkerSize', 10)
% plot(h_vec_km,error_guess,'-o');
% legend('Simulated','Delta-V guess')
% xlabel('Target orbit altitude (km)')
% ylabel('error (%)')
% plot(h_vec_km(bad),error(bad),'rx','MarkerSize', 10)
% % yline(0,'g');
% % yline(10,'r');
% % yline(-10,'r');
% hold off

 


% figure 
% hold on; grid on
% for k = 3:3
%     if exitflag_vec(k) > 0
%     plot(data_all(k).t,(data_all(k).h/1000),'DisplayName', sprintf('h=%.0f km', h_vec(k)/1000));
%     end
% end
% xlabel('t (s)')
% ylabel('altitude (km)')
% title('Electron mission profile (orbit altitude=600 km)')
% % legend('Location','bestoutside')
% hold off

%% Script functions

function [c,ceq] = constraints(sol,sol_scale,h_orbit,lat_launch,azimuth)

global mu r_earth 

sol_unscaled=sol.*sol_scale';

data =sim(sol_unscaled,lat_launch,azimuth);

[h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit, sol_unscaled);

h_min = 100e3;
h_cut_min = 150e3;         
h_SECO=data.h(end);
e_max=0.999;

x = data.x(end);
z = data.z(end);
vx = data.vx(end);
vz = data.vz(end);

r  = hypot(x,z);
urx = x/r;  
urz = z/r;
utx = urz;
utz = -urx;
vr = vx*urx + vz*urz;
vt = vx*utx + vz*utz;
theta_max = deg2rad(5);      
k = tan(theta_max);

v_ref = sqrt(mu/(r_earth + h_orbit));   % ~orbital speed at target radius

theta_constraints = [
    (vr - k*vt)/v_ref;       % vr <= k*vt
   (-vr - k*vt)/v_ref       % -vr <= k*vt  -> |vr| <= k*vt
];

    % equality constraints: hit target altitude, mass
    ceq =  [(h_apo - h_orbit)/h_orbit;
           ];

    c = [(m_prop_2b-m_prop_remaining)/100;
         (h_min-h_per)/h_orbit;
         (h_cut_min - h_SECO)/h_orbit;
         e - e_max;
         theta_constraints
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
alpha_0=sol(6);
alpha_rate=sol(7);

opts_ode = odeset('RelTol',1e-6,'AbsTol',1e-6);
   
% Burn 1a (pitchover)

m0=stage(1).m_dry + stage(1).m_prop + stage(2).m_dry + stage(2).m_prop + m_pld;
m_dot_1=-stage(1).thrust/(stage(1).isp * g);
t_burn_1= -stage(1).m_prop / m_dot_1;
t_burn_1a=t_pitchover+t_pitch;
vx0 = omega_earth * r_earth * cos(lat_launch) * sin(azimuth);

Y0_1a=[0; r_earth; vx0; 0; m0]; %x=0, z=r_earth, vx=depends on location, vz=0, m=GLOW

[t1a,Y1a]=ode113(@(t,Y1a) burn_1a(t,Y1a,stage,m_dot_1,t_pitchover,theta_kick,t_pitch,lat_launch,azimuth),[0 t_burn_1a],Y0_1a,opts_ode);

% Burn 1b (gravity turn)

t_burn_1b=t_burn_1-t_burn_1a;
Y0_1b=Y1a(end,:)'; 

[t1b,Y1b]=ode113(@(t,Y1b) burn_1b(t,Y1b,stage,m_dot_1,lat_launch,azimuth),[t_burn_1a t_burn_1b+t_burn_1a],Y0_1b,opts_ode);

x_1b=Y1b(:,1);
z_1b=Y1b(:,2);
vx_1b=Y1b(:,3);
vz_1b=Y1b(:,4);
m_1b=Y1b(:,5);

x_MECO = x_1b(end);
z_MECO = z_1b(end);
vx_MECO = vx_1b(end);
vz_MECO = vz_1b(end);
m_MECO = m_1b(end) - stage(1).m_dry;

% stage seperation

t_ignition_2=t_burn_1+t_seperation;

Y0_sep=[x_MECO;z_MECO;vx_MECO;vz_MECO;m_MECO]; %state at MECO 

[t_sep,Y_sep]=ode113(@(t,Y_sep) seperation(t,Y_sep,stage,lat_launch,azimuth),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

% Burn 2 (up to apogee) 

m_dot_2= -stage(2).thrust/(stage(2).isp*g);
t_SECO=t_burn_2+t_burn_1+t_seperation;

Y0_2=Y_sep(end,:)';  %state at SES

[t2,Y2]=ode113(@(t,Y2) burn_2(t,Y2,stage,m_dot_2,lat_launch,azimuth,alpha_0,alpha_rate,t_ignition_2),[t_ignition_2 t_SECO],Y0_2,opts_ode);

% Data output

data.t=     [t1a;       t1b(2:end);      t_sep(2:end);      t2(2:end)];
data.x=     [Y1a(:,1);  Y1b(2:end,1);    Y_sep(2:end,1);    Y2(2:end,1)];
data.z=     [Y1a(:,2);  Y1b(2:end,2);    Y_sep(2:end,2);    Y2(2:end,2)];
data.vx=    [Y1a(:,3);  Y1b(2:end,3);    Y_sep(2:end,3);    Y2(2:end,3)];
data.vz=    [Y1a(:,4);  Y1b(2:end,4);    Y_sep(2:end,4);    Y2(2:end,4)];
data.m =    [Y1a(:,5);  Y1b(2:end,5);    Y_sep(2:end,5);    Y2(2:end,5)];

    r = sqrt(data.x.^2 + data.z.^2); 

data.h=r-r_earth;     
    
    v_atm = omega_earth .* r .* cos(lat_launch) .* sin(azimuth);
    urx = data.x./r;  
    urz = data.z./r;
    utx = urz;   
    utz = -urx;
   
    vr = data.vx.*urx + data.vz.*urz;
    vt = data.vx.*utx + data.vz.*utz;

    % idx = vt < 0;
    % vr(idx) = -vr(idx);
    % vt(idx) = -vt(idx);

data.theta= atan2(vr,vt);      
%data.theta = unwrap(data.theta);
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

s = 1e-2; % smoothing scale (m/s) - tune
delta_v_2b = 0.5*((vcirc - va) + sqrt((vcirc - va)^2 + s^2));

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

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

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
    tx =(vx-v_atm)/v_rel; 
    tz =vz/v_rel;
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

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h);
Cd=stage(1).Cd;
A=stage(1).A;

g_x=-(mu*x/r^3);
g_z=-(mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

tx =(vx-v_atm)/v_rel; 
tz =vz/v_rel;
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

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

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

function dY2dt= burn_2(t,Y2,stage,m_dot_2,lat_launch,azimuth,alpha_0,alpha_rate,t_ignition_2)
x=Y2(1);
z=Y2(2);
vx=Y2(3);
vz=Y2(4);
m=Y2(5);

global mu r_earth omega_earth 

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h);
Cd=stage(2).Cd;
A=stage(2).A;

g_x=-(mu*x/r^3);
g_z=-(mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

ux = (vx-v_atm)/v_rel; 
uz = vz/v_rel;
nx = -uz;  
nz = ux;

alpha=alpha_0+alpha_rate*(t-t_ignition_2);

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
