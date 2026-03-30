
close all
clc
clear

vehicle_table=readtable("vehicles_v1.csv");

%% ----Constants----
global g G m_earth r_earth mu h_orbit v_earth stage t_seperation omega_earth

g=9.80665; %units=m/s^2
G=6.674e-11; %units=Nm^2/kg^2
m_earth=5.972e24; %units = kg
r_earth=6371e3; %units = m
mu=G*m_earth;
omega_earth=7.2921159e-5; %units=rads/s

%% ----Atmospheric Model----
h_vac = 84e3;          % model transitions to vacuum above this
dh    = 100;           % grid spacing (m)
h_tab = (0:dh:h_vac)'; % column vector
rho_tab = zeros(size(h_tab));
P_tab   = zeros(size(h_tab));

for i = 1:numel(h_tab)
    [rho_tab(i), P_tab(i)] = isa_model_setup(h_tab(i));  
end

isa.h = h_tab;
isa.rho = rho_tab;
isa.P = P_tab;

%% ----Mission Setup----
h_orbit=400e3; %m
lat_launch=deg2rad(28.5); 
i_target=deg2rad(28.5); %must be greater than abs(lat_launch)

if i_target<(pi/2)
azimuth=asin(cos(i_target)/cos(lat_launch));
else 
    azimuth = pi - asin(cos(i_target)/cos(lat_launch));
end

%% ----Loop----

N = 1;                          %number of rows to simulate?
m_pld_max_vec = nan(N,1);       % payload result (kg)
h_apo_vec     = nan(N,1);       % achieved apogee (m)
h_per_vec     = nan(N,1);       % achieved perigee (m)
iter_time = nan(N,1);

idx = randperm(height(vehicle_table), N);   % unique random rows

for i = 1:N
    
%% Setup

t_iter=tic;
k=idx(i);
stage=vehicle_table(k,:);

v_earth=omega_earth*r_earth*cos(lat_launch); %pure eastward velocity
v_orbit=sqrt(G*m_earth/(r_earth+h_orbit));

t_seperation=10;
burn_time_max=stage.m_prop_2 /(stage.T_stage2/(stage.isp_2*g));
burn_time_guess=300+((burn_time_max-300)/2);

pld_guess=stage.m_wet*0.01;


%% Optimisation

%scale variables to around 1 
u_scale = [deg2rad(5); 10; burn_time_guess; 10; deg2rad(6); deg2rad(0.05)];

%initial guesses
u_0=   [ 
    deg2rad(1); %theta_kick
    10; %t_pitchover
    burn_time_guess; %burn_time_2
    10; %t_pitch
    deg2rad(0); %alpha_0
    deg2rad(0) %alpha rate
];

lowbound=[ 
    deg2rad(0.5); %theta_kick
    5; %t_pitchover
    250; %burn_time_2
    5; %t_pitch
    deg2rad(-6); %alpha_0
    deg2rad(-0.05) %alpha rate
];

upbound= [ 
    deg2rad(5); %theta_kick
    30;  %t_pitchover
    burn_time_max;  %burn_time_2
    20;  %t_pitch
    deg2rad(6); %alpha_0
    deg2rad(0.05) %alpha rate
]; 


u_options=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'FiniteDifferenceType','central' , ...
    'Display','iter', ...
    'FiniteDifferenceStepSize',1e-3, ...
    'ConstraintTolerance',2e-1, ...
    'OptimalityTolerance', 1e-3, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5e5 ...
    );

tol_ceq_km = 0.25;   % apoapsis error allowed
tol_c   = 2e-3;   % small slack on inequalities

m_start = 0.60*pld_guess;     % conservative start
growth  = 1.10;               % +10% per step
maxGrow = 30;

% seed solution (payload + guidance)
sol_seed = [m_start; u_0];

[sol_best, feas] = solve_guidance(m_start, sol_seed, u_scale, lowbound, upbound, ...
    h_orbit, lat_launch, azimuth, isa, tol_ceq_km, tol_c, u_options);

if ~feas
    % If this happens, start even lower (e.g. 0.1*pld_guess) or widen bounds.
    error("Feasibility-first could not find a feasible trajectory at low payload.");
end

m_lo = m_start;
m_hi = NaN;

% Grow payload until it fails (find bracket)
for kk = 1:maxGrow
    m_try = m_lo * growth;

    [sol_try, feas] = solve_guidance( ...
        m_try, sol_best, u_scale, lowbound, upbound, ...
        h_orbit, lat_launch, azimuth, isa, tol_ceq_km, tol_c, u_options);

    if feas
        m_lo = m_try;
        sol_best = sol_try;      % warm start continuation
    else
        m_hi = m_try;
        break
    end
end

% If we never failed, accept m_lo as "at least this much"
if isnan(m_hi)
    sol_optimal = sol_best;
else
    % Bisection to find max feasible payload
    tol_payload = 1.0;   % kg resolution
    maxBisect   = 25;

    for it = 1:maxBisect
        m_mid = 0.5*(m_lo + m_hi);

        [sol_mid, feas] = solve_guidance( ...
            m_mid, sol_best, u_scale, lowbound, upbound, ...
            h_orbit, lat_launch, azimuth, isa, tol_ceq_km, tol_c, u_options);

        if feas
            m_lo = m_mid;
            sol_best = sol_mid;  % keep best feasible solution for warm starts
        else
            m_hi = m_mid;
        end

        if (m_hi - m_lo) < tol_payload
            break
        end
    end

    sol_optimal = sol_best;
end


%% Results 

data=sim(sol_optimal,lat_launch,azimuth,isa);
[h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit,sol_optimal);

m_reserve=(m_prop_remaining - m_prop_2b);
m_pld_max = sol_optimal(1);

% store
vehicle_table.max_payload(k) = m_pld_max;
h_apo_vec(i)     = h_apo;
h_per_vec(i)     = h_per;
iter_time(i) = toc(t_iter);  % seconds elapsed
fprintf("k=%d/%d: %.2f s\n", i, N, iter_time(i));

end


%% ----Script functions----


function [sol_out, feas] = solve_guidance( m_pld, sol_start, u_scale, lowbound, upbound, ...
    h_orbit, lat_launch, azimuth, isa, tol_ceq_km, tol_c, u_options)

    % Extract guidance start
    u_0  = sol_start(2:7);

    % Scale guidance vars (elementwise)
    u_0_s = u_0 ./ u_scale;
    lb_s = lowbound ./ u_scale;
    ub_s = upbound ./ u_scale;

    % Feasibility penalty objective (guidance-only)
    objA = @(us) u_feasibility(us, u_scale, m_pld, h_orbit, lat_launch, azimuth, isa);

    try
        [u_star_s, ~, exitflag] = fmincon(objA, u_0_s, [],[],[],[], lb_s, ub_s, [], u_options);
        u_star = u_star_s .* u_scale;
        sol_out = [m_pld; u_star];

        % Hard feasibility check against your true constraints
        [c, ceq] = constraints_unscaled(sol_out, h_orbit, lat_launch, azimuth, isa);

        feas = (exitflag > 0) && all(isfinite([c; ceq])) ...
            && (abs(ceq(1)) <= tol_ceq_km) && (max(c) <= tol_c);

    catch
        sol_out = sol_start;  % return something safe
      
        feas = false;
    end
end

function [c,ceq] = constraints_unscaled(sol, h_orbit, lat_launch, azimuth, isa)

global mu r_earth 

data =sim(sol,lat_launch,azimuth,isa);

[h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit, sol);

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
    ceq =  [(h_apo - h_orbit)/1000;
           ];

    c = [(m_prop_2b-m_prop_remaining)/10;
         (h_min-h_per)/1000;
         (h_cut_min - h_SECO)/1000;
         e - e_max;
         theta_constraints
    ];
end 

function[rho_atm,P_atm] = isa_model_setup(h)

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

function [rho_atm,P_atm] = isa_model(h, isa)
   
if h >= isa.h(end)
        rho_atm = 0;
        P_atm   = 0;
        return
end
    hq = max(h,0);
    rho_atm = interp1(isa.h, isa.rho, hq, 'linear');
    P_atm   = interp1(isa.h, isa.P,   hq, 'linear');

end

function data =sim(sol,lat_launch,azimuth,isa)

global g stage t_seperation r_earth omega_earth

m_pld=sol(1);
theta_kick=sol(2);
t_pitchover=sol(3);
t_burn_2=sol(4);
t_pitch=sol(5);
alpha_0=sol(6);
alpha_rate=sol(7);

opts_ode = odeset('RelTol',1e-4,'AbsTol',1e-4);
   
% Burn 1a (pitchover)

m0=stage.m_wet + m_pld;
m_dot_1=-stage.T_stage1/(stage.isp_1 * g);
t_burn_1= -stage.m_prop_1/m_dot_1;
t_burn_1a=t_pitchover+t_pitch;
vx0 = omega_earth * r_earth * cos(lat_launch) * sin(azimuth);

Y0_1a=[0; r_earth; vx0; 0; m0]; %x=0, z=r_earth, vx=depends on location, vz=0, m=GLOW

[t1a,Y1a]=ode113(@(t,Y1a) burn_1a(t,Y1a,m_dot_1,t_pitchover,theta_kick,t_pitch,lat_launch,azimuth,isa),[0 t_burn_1a],Y0_1a,opts_ode);

% Burn 1b (gravity turn)

t_burn_1b=t_burn_1-t_burn_1a;
Y0_1b=Y1a(end,:)'; 

[t1b,Y1b]=ode113(@(t,Y1b) burn_1b(t,Y1b,m_dot_1,lat_launch,azimuth,isa),[t_burn_1a t_burn_1b+t_burn_1a],Y0_1b,opts_ode);

x_1b=Y1b(:,1);
z_1b=Y1b(:,2);
vx_1b=Y1b(:,3);
vz_1b=Y1b(:,4);
m_1b=Y1b(:,5);

x_MECO = x_1b(end);
z_MECO = z_1b(end);
vx_MECO = vx_1b(end);
vz_MECO = vz_1b(end);
m_MECO = m_1b(end) - stage.m_dry_1;

% stage seperation

t_ignition_2=t_burn_1+t_seperation;

Y0_sep=[x_MECO;z_MECO;vx_MECO;vz_MECO;m_MECO]; %state at MECO 

[t_sep,Y_sep]=ode113(@(t,Y_sep) seperation(t,Y_sep,lat_launch,azimuth,isa),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

% Burn 2 (up to apogee) 

m_dot_2= -stage.T_stage2/(stage.isp_2*g);
t_SECO=t_burn_2+t_burn_1+t_seperation;

Y0_2=Y_sep(end,:)';  %state at SES

[t2,Y2]=ode113(@(t,Y2) burn_2(t,Y2,m_dot_2,lat_launch,azimuth,alpha_0,alpha_rate,t_ignition_2,isa),[t_ignition_2 t_SECO],Y0_2,opts_ode);

% Data output

data.t=     [t1a;       t1b(2:end);      t_sep(2:end);      t2(2:end)];
data.x=     [Y1a(:,1);  Y1b(2:end,1);    Y_sep(2:end,1);    Y2(2:end,1)];
data.z=     [Y1a(:,2);  Y1b(2:end,2);    Y_sep(2:end,2);    Y2(2:end,2)];
data.vx=    [Y1a(:,3);  Y1b(2:end,3);    Y_sep(2:end,3);    Y2(2:end,3)];
data.vz=    [Y1a(:,4);  Y1b(2:end,4);    Y_sep(2:end,4);    Y2(2:end,4)];
data.m =    [Y1a(:,5);  Y1b(2:end,5);    Y_sep(2:end,5);    Y2(2:end,5)];
data.h=     sqrt(data.x.^2 + data.z.^2)-r_earth;     
    
    r = sqrt(data.x.^2 + data.z.^2); 
    v_atm = omega_earth .* r .* cos(lat_launch) .* sin(azimuth);
    urx = data.x./r;  
    urz = data.z./r;
    utx = urz;   
    utz = -urx;
   
    vr = data.vx.*urx + data.vz.*urz;
    vt = data.vx.*utx + data.vz.*utz;

data.theta= atan2(vr,vt);      
data.v_inf= sqrt((data.vx-v_atm).^2 + data.vz.^2);

end 

function [h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, ~ ,sol)

 global mu r_earth g stage

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

m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage.isp_2*g)));
m_prop_remaining=m_SECO-m_pld-stage.m_dry_2;

end

function J = u_feasibility(u_scaled, u_scale, m_pld, h_orbit, lat_launch, azimuth, isa)

    u = (u_scaled .* u_scale);  % unscale
    sol = [m_pld; u];

    [c, ceq] = constraints_unscaled(sol, h_orbit, lat_launch, azimuth, isa);

    % squared apoapsis error
    J = 1e3*(ceq(1)^2);

    % hinge penalty for inequalities: only penalize when violated (c>0)
    cv = max(c,0);
    J = J + 1e3*sum(cv.^2);

end 

%% ----ODE functions----

function dY1adt= burn_1a(t,Y1a,m_dot_1,t_pitchover,theta_kick,t_pitch,lat_launch,azimuth,isa)
x=Y1a(1);
z=Y1a(2);
vx=Y1a(3);
vz=Y1a(4);
m=Y1a(5);

global mu r_earth omega_earth stage 

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.3;
A=stage.A_1;

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

T=stage.T_stage1-(P_atm*(stage.A_exit_1));
T_x=T*tx;
T_z=T*tz;


dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_1;

dY1adt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY1bdt= burn_1b(~,Y1b,m_dot_1,lat_launch,azimuth,isa)

x=Y1b(1);
z=Y1b(2);
vx=Y1b(3);
vz=Y1b(4);
m=Y1b(5);

global mu r_earth omega_earth stage

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.3;
A=stage.A_1;

g_x=-(mu*x/r^3);
g_z=-(mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

tx =(vx-v_atm)/v_rel; 
tz =vz/v_rel;
T=stage.T_stage1-(P_atm*(stage.A_exit_1));
T_x=T*tx;
T_z=T*tz;

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_1;

dY1bdt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY_sepdt= seperation(~,Y_sep,lat_launch,azimuth,isa)
x=Y_sep(1);
z=Y_sep(2);
vx=Y_sep(3);
vz=Y_sep(4);
m=Y_sep(5);

global mu r_earth omega_earth stage

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,~]=isa_model(h,isa);
Cd=0.3;
A=stage.A_2;

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

function dY2dt= burn_2(t,Y2,m_dot_2,lat_launch,azimuth,alpha_0,alpha_rate,t_ignition_2,isa)
x=Y2(1);
z=Y2(2);
vx=Y2(3);
vz=Y2(4);
m=Y2(5);

global mu r_earth omega_earth stage

r=sqrt(x^2 + z^2);
h=r-r_earth;

v_atm = omega_earth * r * cos(lat_launch) * sin(azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.3;
A=stage.A_2;

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
T=stage.T_stage2-(P_atm*(stage.A_exit_2));
T_x=T*tx;
T_z=T*tz;

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_2;

dY2dt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

    end
