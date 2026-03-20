
close all
clc
clear

vehicle_table=readtable("vehicles_v5_drl_alu.csv");

%% ----Constants----

const.g=9.80665; %units=m/s^2
const.G=6.674e-11; %units=Nm^2/kg^2
const.m_earth=5.972e24; %units = kg
const.r_earth=6371e3; %units = m
const.mu=const.G*const.m_earth;
const.omega_earth=7.2921159e-5; %units=rads/s
const.t_seperation=10;

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

F_rho = griddedInterpolant(isa.h, isa.rho, 'linear', 'nearest');
F_P   = griddedInterpolant(isa.h, isa.P,   'linear', 'nearest');

isa.F_rho = F_rho;
isa.F_P   = F_P;

%% ----Mission Setup----
h_orbit=400e3; %m
const.lat_launch=deg2rad(5); 
i_target=deg2rad(5); %must be greater than abs(lat_launch)

if i_target<(pi/2)
const.azimuth=asin(cos(i_target)/cos(const.lat_launch));
else 
    const.azimuth = pi - asin(cos(i_target)/cos(const.lat_launch));
end

dV_orbit=sqrt(const.mu/(const.r_earth+h_orbit)); %m/s 
dV_grav_loss=0.6*sqrt((2*const.g*h_orbit)/(1+(h_orbit/const.r_earth))); %m/s (approx.)
dV_drag=400; %m/s (approx.)
v_earth=const.omega_earth*const.r_earth*cos(const.lat_launch); %pure eastward velocity
dV_req=dV_orbit+dV_drag+dV_grav_loss-v_earth;

%% ----Loop----

N = height(vehicle_table);                        %number of rows to simulate?
h_apo_vec     = nan(N,1);       % achieved apogee (m)
exitflag_vec = nan(N,1);
iter_time = nan(N,1);
pld_vec = nan(N,1);
guess_pld_vec = nan(N,1);

%rng(1, "twister"); 
idx = randperm(height(vehicle_table), N);   % 10 unique random rows

for i = 1:N
    
%% Setup

t_iter=tic;
% k=idx(i);
stage=vehicle_table(i,:);

v_orbit=sqrt(const.G*const.m_earth/(const.r_earth+h_orbit));

burn_time_max=stage.burn_time_2;
burn_time_guess=150+((burn_time_max-150)/2);

m_wet_1=stage.m_dry_1+stage.m_prop_1;
m_wet_2=stage.m_dry_2+stage.m_prop_2;

dv_residual = @(m_pld) ...
    stage.isp_1 * const.g * log((m_wet_1 + m_wet_2 + m_pld) / ...
                           (stage.m_dry_1 + m_wet_2 + m_pld)) + ...
    stage.isp_2 * const.g * log((m_wet_2 + m_pld) / ...
                           (stage.m_dry_2 + m_pld)) ...
    - dV_req ;

m_pld_lower = 1; % kg minimum
m_pld_upper = 0.2 * stage.m_wet;

if dv_residual(m_pld_lower) * dv_residual(m_pld_upper) <= 0
    pld_guess = fzero(dv_residual, [m_pld_lower, m_pld_upper]);
else
    % Fallbacks if the rocket is drastically under/over-powered
    if dv_residual(m_pld_lower) < 0
        pld_guess = 10; % Vehicle cannot reach orbit even empty
    else
        pld_guess = m_pld_upper; % Super heavy lift vehicle
    end
end

sim_cached([], [], [], [],[], true);

%% Optimisation

%scale variables to around 1 
sol_scale = [pld_guess, deg2rad(5), 10, burn_time_guess, 10, deg2rad(6), deg2rad(0.05)];

%initial guesses
sol_0=   [ 
    0.5*pld_guess; %m_pld
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
    150; %burn_time_2
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
objective = @(sol) obj_payload_leftover(sol, sol_scale, lambda, isa,stage,const);

nonlcon=@(sol) constraints(sol,sol_scale,h_orbit,isa,stage,const);

iter_start = tic;

options=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'FiniteDifferenceType','forward' , ...
    'Display','final', ...
    'FiniteDifferenceStepSize',1e-4, ...
    'ConstraintTolerance',5e-2, ...
    'OptimalityTolerance', 1e-3, ...
    'MaxIterations', 50, ...
    'MaxFunctionEvaluations', 2000 ...
    );

[sol__optimal_scaled, fval, exitflag, output]=fmincon(objective,sol_0,[],[],[],[],lowbound,upbound,nonlcon,options);
sol_optimal=sol__optimal_scaled.*sol_scale';

%% Results 

data=sim(sol_optimal,isa,stage,const);
[h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, h_orbit,sol_optimal,stage,const);

m_reserve=(m_prop_remaining - m_prop_2b);
m_pld_max = sol_optimal(1);
error=(h_apo-h_orbit)/h_orbit;

% store
if exitflag > 0 && abs(error) < 0.1
vehicle_table.max_payload(i) = m_pld_max;

else 
    vehicle_table.max_payload(i) = 1;
end

pld_vec(i)       = (m_pld_max/stage.m_wet);
guess_pld_vec(i)       = (pld_guess/stage.m_wet);
h_apo_vec(i)     = h_apo;
exitflag_vec(i) = exitflag;
iter_time(i) = toc(t_iter);  % seconds elapsed
fprintf("k=%d/%d: %.2f s\n", i, N, iter_time(i));

end

av_iter_time=mean(iter_time);
ok = exitflag_vec > 0;
rate_success=numel(exitflag_vec(ok))./N;
fprintf('av_iter_time = %.2f s\n', av_iter_time)
fprintf('av_exitflag %.2f %\n', rate_success*100)

%% --Write CSV--

outFile = "Kourou_5deg_400km_rtls_alu.csv";
writetable(vehicle_table, outFile);

fprintf("Wrote: %s\n", outFile);

%% ----Script functions----

function [c,ceq] = constraints(sol,sol_scale,h_orbit,isa,stage,const)

    try
    cache = sim_cached(sol, sol_scale, isa, stage, const, false);

    % If anything important is non-finite, force infeasible
    if ~all(isfinite([cache.h_apo, cache.h_per, cache.e, cache.m_prop_2b, cache.m_prop_remaining])) ...
            || ~all(isfinite([cache.data.x(end), cache.data.z(end), cache.data.vx(end), cache.data.vz(end)]))
        c   = 1e6 * ones(6,1);   % positive => violated
        ceq = 1e6;              % not equal to 0 => violated
        return
    end

data = cache.data;
h_apo = cache.h_apo;
h_per = cache.h_per;
e = cache.e;
m_prop_2b = cache.m_prop_2b;
m_prop_remaining = cache.m_prop_remaining;

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

v_ref = sqrt(const.mu/(const.r_earth + h_orbit));   % ~orbital speed at target radius

theta_constraints = [
   (vr - k*vt)/v_ref;       % vr <= k*vt
   (-vr - k*vt)/v_ref       % -vr <= k*vt  -> |vr| <= k*vt
];

    % equality constraints: hit target altitude, mass
    ceq =  [(h_apo - h_orbit)/h_orbit;
           ];

    c = [(m_prop_2b-m_prop_remaining)/stage.m_prop_2;
         (h_min-h_per)/h_orbit;
         (h_cut_min - h_SECO)/h_orbit;
         e - e_max;
         theta_constraints
    ];
    
    catch
  c   = NaN(6,1);
        ceq = NaN;
    end
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
    rho_atm = isa.F_rho(hq);
    P_atm   = isa.F_P(hq);

end

function data =sim(sol,isa,stage,const)

m_pld=sol(1);
theta_kick=sol(2);
t_pitchover=sol(3);
t_burn_2=sol(4);
t_pitch=sol(5);
alpha_0=sol(6);
alpha_rate=sol(7);

opts_ode = odeset('RelTol',1e-5,'AbsTol',1e-5,'MaxStep',50);  

% Burn 1a (pitchover)
m0=stage.m_wet + m_pld;
m_dot_1=-stage.T_stage1/(stage.isp_1 * const.g);
t_burn_1= -(stage.m_prop_1*0.85)/m_dot_1;
t_burn_1a=t_pitchover+t_pitch;
vx0 = const.omega_earth * const.r_earth * cos(const.lat_launch) * sin(const.azimuth);

Y0_1a=[0; const.r_earth; vx0; 0; m0]; %x=0, z=r_earth, vx=depends on location, vz=0, m=GLOW

[~,Y1a]=ode113(@(t,Y1a) burn_1a(t,Y1a,m_dot_1,t_pitchover,theta_kick,t_pitch,isa,stage,const),[0 t_burn_1a],Y0_1a,opts_ode);

% Burn 1b (gravity turn)
t_burn_1b=t_burn_1-t_burn_1a;
Y0_1b=Y1a(end,:)'; 

[~,Y1b]=ode113(@(t,Y1b) burn_1b(t,Y1b,m_dot_1,isa,stage,const),[t_burn_1a t_burn_1b+t_burn_1a],Y0_1b,opts_ode);

x_1b=Y1b(:,1);
z_1b=Y1b(:,2);
vx_1b=Y1b(:,3);
vz_1b=Y1b(:,4);
m_1b=Y1b(:,5);

x_MECO = x_1b(end);
z_MECO = z_1b(end);
vx_MECO = vx_1b(end);
vz_MECO = vz_1b(end);
m_MECO = stage.m_dry_2 + stage.m_prop_2 + m_pld;

% stage seperation
t_ignition_2=t_burn_1+const.t_seperation;

Y0_sep=[x_MECO;z_MECO;vx_MECO;vz_MECO;m_MECO]; %state at MECO 

[~,Y_sep]=ode113(@(t,Y_sep) seperation(t,Y_sep,isa,stage,const),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

% Burn 2 (up to apogee) 
m_dot_2= -stage.T_stage2/(stage.isp_2*const.g);
t_SECO=t_burn_2+t_burn_1+const.t_seperation;

Y0_2=Y_sep(end,:)';  %state at SES

[t2,Y2]=ode113(@(t,Y2) burn_2(t,Y2,m_dot_2,alpha_0,alpha_rate,t_ignition_2,isa,stage,const),[t_ignition_2 t_SECO],Y0_2,opts_ode);

% Data output
data.t=     [t2(end)];
data.x=     [Y2(end,1)];
data.z=     [Y2(end,2)];
data.vx=    [Y2(end,3)];
data.vz=    [Y2(end,4)];
data.m =    [Y2(end,5)];
data.h=     sqrt(data.x.^2 + data.z.^2)-const.r_earth;     

end 

function [h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, ~ ,sol,stage,const)

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

orbit_energy_sp  = v0mag^2/2 - const.mu/r0mag;

if ~isfinite(orbit_energy_sp) || abs(orbit_energy_sp) < 1e-12
    error("Bad orbit energy");
end

h_vec = r0(1)*v0(2) - r0(2)*v0(1);
h_abs = abs(h_vec);

e = sqrt(1 + 2*orbit_energy_sp*h_abs^2/const.mu^2);
a = -const.mu/(2*orbit_energy_sp);

r_per = a*(1 - e);
r_apo = a*(1 + e);

h_per = r_per - const.r_earth;
h_apo = r_apo - const.r_earth;

% circularization at apogee 
va = sqrt(const.mu*(2/r_apo - 1/a));
vcirc = sqrt(const.mu/r_apo);

s = 0.5; % smoothing scale (m/s) - tune
delta_v_2b = 0.5*((vcirc - va) + sqrt((vcirc - va)^2 + s^2));

m_prop_2b=m_SECO*(1-exp(-delta_v_2b/(stage.isp_2*const.g)));
m_prop_remaining=m_SECO-m_pld-stage.m_dry_2;

end

function J = obj_payload_leftover(sol, sol_scale, lambda, isa,stage,const)

    try
cache=sim_cached(sol, sol_scale, isa,stage,const,false);

    % If anything is non-finite, penalize
    if ~isfinite(cache.m_prop_remaining) || ~isfinite(cache.m_prop_2b) || ~isfinite(cache.sol(1))
        J = 1e12;
        return
    end

m_pld  = cache.sol(1);     
x = cache.m_prop_remaining - cache.m_prop_2b;  % leftover after planned circ (kg)
s = 1e-2;                 % kg scale, tune
m_left = 0.5*(x + sqrt(x.^2 + s^2)); 

% maximize (m_pld - lambda*m_left)  <=>  minimize negative:
J = -(m_pld - lambda*m_left);

    if ~isfinite(J)
        J = 1e12;
    end

    catch
    J = 1e12;  % bad point -> huge penalty
     end
end

function cache = sim_cached(sol, sol_scale, isa, stage, const, reset)

    persistent last_sol_scaled last_out

 if nargin >= 6 && reset
        last_sol_scaled = [];
        last_out = [];
        cache = [];
        return
 end

    % First call or new design point?
   tol = 1e-8;
    if isempty(last_sol_scaled) || norm(sol - last_sol_scaled) > tol
        sol_unscaled = sol(:) .* sol_scale(:);

        data = sim(sol_unscaled, isa, stage,const);

if ~all(isfinite([data.x(end), data.z(end), data.vx(end), data.vz(end), data.m(end), data.h(end)]))
    error("sim returned non-finite state");
end

        % Precompute orbit quantities once too (so you don't redo in both places)
        [h_apo, h_per, e, m_prop_2b, m_prop_remaining] = orbit_calc(data, [], sol_unscaled,stage,const);

        % Package everything objective/constraints need
        cache.sol = sol_unscaled;
        cache.data = data;
        cache.h_apo = h_apo;
        cache.h_per = h_per;
        cache.e = e;
        cache.m_prop_2b = m_prop_2b;
        cache.m_prop_remaining = m_prop_remaining;

        % Save to cache
        last_sol_scaled = sol;
        last_out = cache;
    else
        cache = last_out;
    end
end

%% ----ODE functions----

function dY1adt= burn_1a(t,Y1a,m_dot_1,t_pitchover,theta_kick,t_pitch,isa,stage,const)
x=Y1a(1);
z=Y1a(2);
vx=Y1a(3);
vz=Y1a(4);
m=Y1a(5);

r=sqrt(x^2 + z^2);
h=r-const.r_earth;

v_atm = const.omega_earth * r * cos(const.lat_launch) * sin(const.azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.2;
A=stage.A_1;

g_x=-(const.mu*x/r^3);
g_z=-(const.mu*z/r^3);

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

T=(stage.T1-(P_atm*(stage.A_exit_1)))*stage.N1;
T_x=T*tx;
T_z=T*tz;

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_1;

dY1adt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY1bdt= burn_1b(~,Y1b,m_dot_1,isa,stage,const)

x=Y1b(1);
z=Y1b(2);
vx=Y1b(3);
vz=Y1b(4);
m=Y1b(5);

r=sqrt(x^2 + z^2);
h=r-const.r_earth;

v_atm = const.omega_earth * r * cos(const.lat_launch) * sin(const.azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.2;
A=stage.A_1;

g_x=-(const.mu*x/r^3);
g_z=-(const.mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

tx =(vx-v_atm)/v_rel; 
tz =vz/v_rel;
T=(stage.T1-(P_atm*(stage.A_exit_1)))*stage.N1;
T_x=T*tx;
T_z=T*tz;

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_1;

dY1bdt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY_sepdt= seperation(~,Y_sep,isa,stage,const)
x=Y_sep(1);
z=Y_sep(2);
vx=Y_sep(3);
vz=Y_sep(4);
m=Y_sep(5);

r=sqrt(x^2 + z^2);
h=r-const.r_earth;

v_atm = const.omega_earth * r * cos(const.lat_launch) * sin(const.azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,~]=isa_model(h,isa);
Cd=0.2;
A=stage.A_2;

g_x=-(const.mu*x/r^3);
g_z=-(const.mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(-D_x)/m;
dvzdt=g_z+(-D_z)/m;
dmdt=0;

dY_sepdt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

end

function dY2dt= burn_2(t,Y2,m_dot_2,alpha_0,alpha_rate,t_ignition_2,isa,stage,const)
x=Y2(1);
z=Y2(2);
vx=Y2(3);
vz=Y2(4);
m=Y2(5);


r=sqrt(x^2 + z^2);
h=r-const.r_earth;

v_atm = const.omega_earth * r * cos(const.lat_launch) * sin(const.azimuth);  % in-plane

v_rel=max(sqrt((vx-v_atm)^2 + vz^2),0.1);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.2;
A=stage.A_2;

g_x=-(const.mu*x/r^3);
g_z=-(const.mu*z/r^3);

D_x=(0.5*rho_atm*v_inf^2*Cd*A)*((vx-v_atm)/v_inf);
D_z=(0.5*rho_atm*v_inf^2*Cd*A)*(vz/v_inf);

ux = (vx-v_atm)/v_rel; 
uz = vz/v_rel;
nx = -uz;  
nz = ux;

alpha=alpha_0+alpha_rate*(t-t_ignition_2);

tx=cos(alpha)*ux + sin(alpha)*nx;
tz=cos(alpha)*uz + sin(alpha)*nz;
T=stage.T_stage2;
T_x=T*tx;
T_z=T*tz;

dxdt=vx;
dzdt=vz;
dvxdt=g_x+(T_x-D_x)/m;
dvzdt=g_z+(T_z-D_z)/m;
dmdt=m_dot_2;

dY2dt = [dxdt;dzdt;dvxdt;dvzdt;dmdt];

    end
