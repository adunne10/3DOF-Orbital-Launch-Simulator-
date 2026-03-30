
close all
clc
clear

vehicle_table=readtable("vehicles_v2.csv");

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
const.lat_launch=deg2rad(28.5); 
i_target=deg2rad(51.6); %must be greater than abs(lat_launch)

if i_target<(pi/2)
const.azimuth=asin(cos(i_target)/cos(const.lat_launch));
else 
    const.azimuth = pi - asin(cos(i_target)/cos(const.lat_launch));
end

%% ----Loop----

N = 30;                        %number of rows to simulate?
h_apo_vec     = nan(N,1);       % achieved apogee (m)
h_per_vec     = nan(N,1);       % achieved perigee (m)
exitflag_vec = nan(N,1);
iter_time = nan(N,1);
pld_vec = nan(N,1);

% rng(3, "twister"); 
idx = randperm(height(vehicle_table), N);   % 10 unique random rows

for i = 1:N
    
%% Setup

t_iter=tic;
k=idx(i);
stage=vehicle_table(k,:);

pld_guess=stage.m_wet*0.02;

%% Optimisation

sim_cached([], [], [], [], [], true);

%scale variables to around 1 
sol_scale = [pld_guess, deg2rad(5), 10, 10, deg2rad(6), deg2rad(0.05)];

%initial guesses
sol_0=   [ 
    pld_guess; %m_pld
    deg2rad(1); %theta_kick
    10; %t_pitchover
    10; %t_pitch
    deg2rad(0); %alpha_0
    deg2rad(0.01) %alpha rate
]./sol_scale';

lowbound=[ 
    (0.5*pld_guess); %m_pld
    deg2rad(0.5); %theta_kick
    5; %t_pitchover
    5; %t_pitch
    deg2rad(-6); %alpha_0
    deg2rad(-0.05) %alpha rate
]./sol_scale';

upbound= [ 
    (10*pld_guess); %m_pld
    deg2rad(5); %theta_kick
    30;  %t_pitchover
    20;  %t_pitch
    deg2rad(6); %alpha_0
    deg2rad(0.05) %alpha rate
]./sol_scale'; 


nonlcon=@(sol) constraints(sol,sol_scale,h_orbit,isa,stage,const);

options=optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'ScaleProblem', 'obj-and-constr', ... 
    'FiniteDifferenceType','forward' , ...
    'Display','final', ...
    'FiniteDifferenceStepSize',1e-2, ...
    'ConstraintTolerance',1e-1, ...
    'OptimalityTolerance', 1e-3, ...
    'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 5e5 ...
    );

objective = @(sol) -(sol(1)*sol_scale(1));   % maximize real payload (kg)
[sol__optimal_scaled, fval, exitflag, output]=fmincon(objective,sol_0,[],[],[],[],lowbound,upbound,nonlcon,options);
sol_optimal=sol__optimal_scaled.*sol_scale';


%% Results 

m_pld_max = sol_optimal(1);

% store
if exitflag > 0 
vehicle_table.max_payload(k) = m_pld_max;
end
pld_vec(i)       = (m_pld_max/stage.m_wet);
exitflag_vec(i) = exitflag;
iter_time(i) = toc(t_iter);  % seconds elapsed
fprintf("k=%d/%d: %.2f s\n", i, N, iter_time(i));

end

av_iter_time=mean(iter_time);
ok = exitflag_vec > 0;
rate_success=numel(exitflag_vec(ok))./N;
fprintf('av_iter_time = %.2f s\n', av_iter_time)
fprintf('av_exitflag %.2f %\n', rate_success*100)

%% ----Script functions----

function [c,ceq] = constraints(sol,sol_scale,h_orbit,isa,stage,const)

    try
        cache = sim_cached(sol, sol_scale, isa, stage, const, false);
        
        % FIXED: Only check states that actually exist in your new data structure
        if ~all(isfinite([cache.data.x(end), cache.data.z(end), cache.data.vx(end), cache.data.vz(end)]))
            c   = NaN(4,1);  
            ceq = NaN;              
            return
        end
        
        data = cache.data;
        v_circ   = sqrt(const.mu / (const.r_earth + h_orbit));
        h_SECO=data.h(end);
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
        
        c_speed = (abs(vt) - 1.02*v_circ)/v_circ;      % vt <= 1.02 vcirc   
        c_speed2 = (0.98*v_circ - abs(vt))/v_circ;     % vt >= 0.98 vcirc
        v_ref = sqrt(const.mu/(const.r_earth + h_orbit));   
        
        theta_constraints = [
           (vr - k*vt)/v_ref;       
           (-vr - k*vt)/v_ref       
        ];
        
        % equality constraints: hit target altitude
        ceq =  (h_SECO - h_orbit)/1000;
               
        c = [theta_constraints;
             c_speed;
             c_speed2];
        
    catch
        % FIXED: Return NaNs so fmincon doesn't break its line-search
        c   = NaN(4,1);
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
t_pitch=sol(4);
alpha_0=sol(5);
alpha_rate=sol(6);

opts_ode = odeset('RelTol',1e-6,'AbsTol',1e-6);
   
% Burn 1a (pitchover)
m0=stage.m_wet + m_pld;
m_dot_1=-stage.T_stage1/(stage.isp_1 * const.g);
t_burn_1= -stage.m_prop_1/m_dot_1;
t_burn_1a=t_pitchover+t_pitch;
vx0 = const.omega_earth * const.r_earth * cos(const.lat_launch) * sin(const.azimuth);

Y0_1a=[0; const.r_earth; vx0; 0; m0]; %x=0, z=r_earth, vx=depends on location, vz=0, m=GLOW

[~,Y1a]=ode113(@(t,Y1a) burn_1a(t,Y1a,m_dot_1,t_pitchover,theta_kick,t_pitch,isa,stage,const),[0 t_burn_1a],Y0_1a,opts_ode);

if isempty(Y1a) || any(~isfinite(Y1a(end,:)))
    error("ODE integration failed or produced non-finite state");
end

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
m_MECO = m_1b(end) - stage.m_dry_1;

% stage seperation
t_ignition_2=t_burn_1+const.t_seperation;

Y0_sep=[x_MECO;z_MECO;vx_MECO;vz_MECO;m_MECO]; %state at MECO 

[~,Y_sep]=ode113(@(t,Y_sep) seperation(t,Y_sep,isa,stage,const),[t_burn_1 t_ignition_2],Y0_sep,opts_ode);

% Burn 2 (up to apogee) 
m_dot_2= -stage.T_stage2/(stage.isp_2*const.g);
t_burn_2=-stage.m_prop_2/m_dot_2;
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

function out = sim_cached(sol_scaled, sol_scale, isa, stage, const, reset)
    persistent keys outs

    if nargin >= 6 && reset
        keys = [];
        outs = {};
        out = [];
        return
    end

    sol_scaled = sol_scaled(:);
    sol_scale  = sol_scale(:);

    % Quantize scaled sol so tiny FP differences still hit the cache
    q   = 1e-8;
    key = round(sol_scaled/q)*q;

    % --- lookup ---
    if ~isempty(keys)
        % keys is (nVars x nStored)
        diffs = max(abs(keys - key), [], 1);
        j = find(diffs == 0, 1);
        if ~isempty(j)
            out = outs{j};
            return
        end
    end

    % --- compute fresh ---
    sol_unscaled = sol_scaled .* sol_scale;
    data = sim(sol_unscaled, isa, stage, const);

    out = struct('sol', sol_unscaled, 'data', data);

    % --- insert ---
    Kmax = 30;

    if isempty(keys)
        % initialize with correct row count
        keys = key;          % (nVars x 1)
        outs = {out};        % cell array of structs
        return
    end

    if size(keys,2) < Kmax
        keys(:, end+1) = key;
        outs{end+1}    = out;
    else
        % simple FIFO: drop oldest (safe + easy)
        keys(:,1) = [];
        outs(1)   = [];
        keys(:, end+1) = key;
        outs{end+1}    = out;
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

v_rel=sqrt((vx-v_atm)^2 + vz^2+0.1^2);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.3;
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

function dY1bdt= burn_1b(~,Y1b,m_dot_1,isa,stage,const)

x=Y1b(1);
z=Y1b(2);
vx=Y1b(3);
vz=Y1b(4);
m=Y1b(5);

r=sqrt(x^2 + z^2);
h=r-const.r_earth;

v_atm = const.omega_earth * r * cos(const.lat_launch) * sin(const.azimuth);  % in-plane

v_rel=sqrt((vx-v_atm)^2 + vz^2+0.1^2);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.3;
A=stage.A_1;

g_x=-(const.mu*x/r^3);
g_z=-(const.mu*z/r^3);

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

function dY_sepdt= seperation(~,Y_sep,isa,stage,const)
x=Y_sep(1);
z=Y_sep(2);
vx=Y_sep(3);
vz=Y_sep(4);
m=Y_sep(5);


r=sqrt(x^2 + z^2);
h=r-const.r_earth;

v_atm = const.omega_earth * r * cos(const.lat_launch) * sin(const.azimuth);  % in-plane

v_rel=sqrt((vx-v_atm)^2 + vz^2+0.1^2);
v_inf=v_rel+0.1;

[rho_atm,~]=isa_model(h,isa);
Cd=0.3;
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

v_rel=sqrt((vx-v_atm)^2 + vz^2+0.1^2);
v_inf=v_rel+0.1;

[rho_atm,P_atm]=isa_model(h,isa);
Cd=0.3;
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
