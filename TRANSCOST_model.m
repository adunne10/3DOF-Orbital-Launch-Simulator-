clear; clc

%% --Design space setup--

vehicle_payload_data=readtable("Andoya_97.5deg_500km_exp_alu.csv");

% Unpack table variables
ID=vehicle_payload_data.ID;

% M1  = vehicle_payload_data.M1;
F1  = vehicle_payload_data.F1 ;
EC1 = vehicle_payload_data.EC1 ;
T1  = vehicle_payload_data.T1;
N1  = vehicle_payload_data.N1;

% M2  = vehicle_payload_data.M2;
F2  = vehicle_payload_data.F2;
EC2 = vehicle_payload_data.EC2;
T2  = vehicle_payload_data.T2;
N2  = vehicle_payload_data.N2;

T_stage1 = vehicle_payload_data.T_stage1;
T_stage2 = vehicle_payload_data.T_stage2;

m_dry_1=vehicle_payload_data.m_dry_1;
m_dry_2=vehicle_payload_data.m_dry_2;
m_prop_1=vehicle_payload_data.m_prop_1;
m_prop_2=vehicle_payload_data.m_prop_2;
max_payload=vehicle_payload_data.max_payload;

isp_1=vehicle_payload_data.isp_1;
isp_2=vehicle_payload_data.isp_2;

m_wet=m_dry_1+m_dry_2+m_prop_1+m_prop_2;
GLOW=m_wet+max_payload;


%% --Engine Properties table--

cycle_CH4=['Expander,closed','Gas generator','Staged, full-flow','Staged, oxidizer-rich'];
index_CH4 = ['A','B','C','D'];
T_W_CH4= [45 113 163 122];
f_1_CH4= [0.7 0.9 1.4 1.2];
O_F_CH4=3.7;
rho_fuel_CH4=422;
C_fuel_CH4=6.5; %£/kg

cycle_LH2=['Expander,closed','Expander, open','Gas generator','Staged, fuel-rich'];
index_LH2 = ['A','B','C','D'];
T_W_LH2= [50 54 58 47];
f_1_LH2= [1.1 1.1 0.8 1];
O_F_LH2=6;
rho_fuel_LH2=70;
C_fuel_LH2=4.5; %£/kg

cycle_RP1=['Combustion tap-off'	'Electric pump'	'Gas generator'	'Staged, oxidizer-rich'];
index_RP1 = ['A','B','C','D'];
T_W_RP1= [86 72.8 103 83];
f_1_RP1= [0.9 0.7 0.4 0.4];
O_F_RP1=2.7;
rho_fuel_RP1=813;
C_fuel_RP1=1.7; %£/kg 

%% --Extra Design Calcs --

idx1 = double(char(EC1)) - double('A') + 1;
idx2 = double(char(EC2)) - double('A') + 1;

O_F_1 = nan(numel(F1),1);
O_F_2 = nan(numel(F2),1);

% --- Stage 1
i = (F1 == "CH4"); O_F_1(i) = O_F_CH4;
i = (F1 == "LH2"); O_F_1(i) = O_F_LH2;
i = (F1 == "RP1"); O_F_1(i) = O_F_RP1;
% --- Stage 2
i = (F2 == "CH4"); O_F_2(i) = O_F_CH4;
i = (F2 == "LH2"); O_F_2(i) = O_F_LH2;
i = (F2 == "RP1"); O_F_2(i) = O_F_RP1;

% Preallocate
T_W_1 = nan(numel(F1),1);
T_W_2 = nan(numel(F2),1);

% --- Stage 1: T/W = f(F1, EC1)
i = (F1 == "CH4"); T_W_1(i) = T_W_CH4(idx1(i));
i = (F1 == "LH2"); T_W_1(i) = T_W_LH2(idx1(i));
i = (F1 == "RP1"); T_W_1(i) = T_W_RP1(idx1(i));
% --- Stage 2: T/W = f(F2, EC2)
i = (F2 == "CH4"); T_W_2(i) = T_W_CH4(idx2(i));
i = (F2 == "LH2"); T_W_2(i) = T_W_LH2(idx2(i));
i = (F2 == "RP1"); T_W_2(i) = T_W_RP1(idx2(i));

m_engine_1= T1./T_W_1/9.81;
m_engine_2= T2./T_W_2/9.81;

% mass_factors=[0.8 1.0 ];% CFRP,ALu
% mass_factor_1 = nan(numel(M1),1);
% mass_factor_2 = nan(numel(M2),1);
% i = (M1=="CFRP");      mass_factor_1(i) =  mass_factors(1);
% i = (M1=="Aluminium"); mass_factor_1(i) =  mass_factors(2);
% i = (M2=="CFRP");      mass_factor_2(i) =  mass_factors(1);
% i = (M2=="Aluminium"); mass_factor_2(i) =  mass_factors(2);

%% --Cost factors (f_n)--

f_0=1.04^2; %2 stage vehicle
f_3=1;      %team with some related experience
p=0.8;      %learning factor
f_6=1;      %on time schedule
f_7=1;      %no parallel contractors
f_8=0.86;   %for Europe

%% --Cost calculations--

% C_work_year=0.25; %million

% Preallocate
f1_1 = nan(numel(F1),1);
f1_2 = nan(numel(F2),1);

% --- Stage 1: Isp = f(F1, EC1)
i = (F1 == "CH4"); f1_1(i) = f_1_CH4(idx1(i));
i = (F1 == "LH2"); f1_1(i) = f_1_LH2(idx1(i));
i = (F1 == "RP1"); f1_1(i) = f_1_RP1(idx1(i));
% --- Stage 2: Isp = f(F2, EC2)
i = (F2 == "CH4"); f1_2(i) = f_1_CH4(idx2(i));
i = (F2 == "LH2"); f1_2(i) = f_1_LH2(idx2(i));
i = (F2 == "RP1"); f1_2(i) = f_1_RP1(idx2(i));

N_Q=100;
f_2_engines=0.026*(log(N_Q))^2;
DD_engine_1=(227.*(m_engine_1).^0.48).*f1_1.*f_2_engines.*f_3;
DD_engine_2=(227.*(m_engine_2).^0.48).*f1_2.*f_2_engines.*f_3;

%expendable stages
DD_stage_1_ex=100.*(m_dry_1-(m_engine_1.*N1)).^0.555*f_3;
DD_stage_2_ex=100.*(m_dry_2-(m_engine_2.*N2)).^0.555*f_3;
%reusable stages
DD_stage_1_re=803.5.*(m_dry_1-(m_engine_1.*N1)).^0.385*f_3;
DD_stage_2_re=803.5.*(m_dry_2-(m_engine_2.*N2)).^0.385*f_3;

N_range = 1:50;       
nN = numel(N_range);
nCfg = numel(ID);

% Preallocate outputs: [config x N_launches]
FU_total_ex_MYr = nan(nCfg,nN);
ground_ops_ex_MYr =nan(nCfg,nN);
MC_MYr=nan(nCfg,nN);

for k = 1:nN
    N_launches = N_range(k);

f4_engine_1=(N_launches*N1).^(log(p)/log(2));
f4_engine_2=(N_launches*N2).^(log(p)/log(2));
f4_stage=(N_launches).^(log(p)/log(2));

FU_engine_cryo_1=f4_engine_1.*5.16.*(N1).*m_engine_1.^0.45;
FU_engine_stor_1=f4_engine_1.*1.9.*(N1).*m_engine_1.^0.535;
FU_engine_cryo_2=f4_engine_2.*5.16.*(N2).*m_engine_2.^0.45;
FU_engine_stor_2=f4_engine_2.*1.9.*(N2).*m_engine_2.^0.535;

% Preallocate
FU_engine_1 = nan(numel(F1),1);
FU_engine_2 = nan(numel(F2),1);

% --- Stage 1
i = (F1 == "CH4"); FU_engine_1(i) = FU_engine_cryo_1(i);
i = (F1 == "LH2"); FU_engine_1(i) = FU_engine_cryo_1(i);
i = (F1 == "RP1"); FU_engine_1(i) = FU_engine_stor_1(i);
% --- Stage 2
i = (F2 == "CH4"); FU_engine_2(i) = FU_engine_cryo_2(i);
i = (F2 == "LH2"); FU_engine_2(i) = FU_engine_cryo_2(i);
i = (F2 == "RP1"); FU_engine_2(i) = FU_engine_stor_2(i);

FU_stage_cryo_1=f4_stage.*1.30.*(m_dry_1-(m_engine_1.*N1)).^0.65;
FU_stage_stor_1= f4_stage.*0.83.*(m_dry_1-(m_engine_1.*N1)).^0.65;
FU_stage_cryo_2=f4_stage.*1.30.*(m_dry_2-(m_engine_2.*N2)).^0.65;
FU_stage_stor_2=f4_stage.*0.83.*(m_dry_2-(m_engine_2.*N2)).^0.65;

% Preallocate
FU_stage_1 = nan(numel(F1),1);
FU_stage_2 = nan(numel(F2),1);

% --- Stage 1
i = (F1 == "CH4"); FU_stage_1(i) = FU_stage_cryo_1(i); 
i = (F1 == "LH2"); FU_stage_1(i) = FU_stage_cryo_1(i);
i = (F1 == "RP1"); FU_stage_1(i) = FU_stage_stor_1(i);
% --- Stage 2
i = (F2 == "CH4"); FU_stage_2(i) = FU_stage_cryo_2(i);
i = (F2 == "LH2"); FU_stage_2(i) = FU_stage_cryo_2(i);
i = (F2 == "RP1"); FU_stage_2(i) = FU_stage_stor_2(i);

% recovery_1=(1.5/N_launches).*(7*N_launches.^0.7+(m_dry_1+m_prop_1).^0.83);
% recovery_2=(1.5/N_launches).*(7*N_launches.^0.7+(m_dry_2+m_prop_2).^0.83);
% 
% refurb_1=0.11.*FU_engine_1+0.01.*FU_stage_1;
% refurb_2=0.11.*FU_engine_2+0.01.*FU_stage_2;

f_v=1.0;

ground_ops_N=8*(m_wet.^0.67).*(5.^(-0.9))*(2^0.7)*0.5*f_v*f4_stage*f_8;
MC_N=20*0.8*5^(-0.65)*f4_stage*f_8;

% Preallocate
C_prop_1 = nan(numel(F1),1);
C_prop_2 = nan(numel(F2),1);

% --- Stage 1
i = (F1 == "CH4"); C_prop_1(i) = C_fuel_CH4;
i = (F1 == "LH2"); C_prop_1(i) = C_fuel_LH2;
i = (F1 == "RP1"); C_prop_1(i) = C_fuel_RP1;
% --- Stage 2
i = (F2 == "CH4"); C_prop_2(i) = C_fuel_CH4;
i = (F2 == "LH2"); C_prop_2(i) = C_fuel_LH2;
i = (F2 == "RP1"); C_prop_2(i) = C_fuel_RP1;


DD_total_ex_MYr=(f_0.*(DD_stage_1_ex+DD_stage_2_ex+DD_engine_1+DD_engine_2).*f_6.*f_7.*f_8); 

FU_total_ex_N=f_0^2.*(FU_stage_1+FU_stage_2+FU_engine_1+FU_engine_2).*f_8; %wrkyr

FU_total_ex_MYr(:,k)=FU_total_ex_N;
ground_ops_ex_MYr(:,k)=ground_ops_N;
MC_MYr(:,k)=MC_N;

end

FU_ave_ex=mean(FU_total_ex_MYr,2); %MYr
GOps_ave_ex=mean(ground_ops_ex_MYr,2); %MYr
MC_ave_ex=mean(MC_MYr,2); %MYr


C_LOX=0.2; %£/kg

C_prop=(C_prop_1.*(m_prop_1.*(1+O_F_1))+C_LOX.*(m_prop_1-m_prop_1.*(1+O_F_1))...
    +C_prop_2.*(m_prop_2.*(1+O_F_2))+C_LOX.*(m_prop_2-m_prop_2.*(1+O_F_2)) )/1e6;

%% --Create table and filter--

valid=max_payload>100;


vehicle_payload_cost_data = table(ID,F1,EC1,T1,N1,T_stage1,m_dry_1,m_prop_1,isp_1, ...
                F2,EC2,T2,N2,T_stage2,m_dry_2,m_prop_2,isp_2, ...
                m_wet,max_payload,GLOW,DD_total_ex_MYr,C_prop, ...
                FU_ave_ex,GOps_ave_ex,MC_ave_ex);

vehicle_payload_cost_data = vehicle_payload_cost_data(valid,:);


fprintf("Total configurations: %d\n", height(vehicle_payload_data));
fprintf("Converged configurations: %d\n", height(vehicle_payload_cost_data));

%% --Write CSV--

outFile = "COSTED_Andoya_97.5deg_500km_exp_alu.csv";
writetable(vehicle_payload_cost_data, outFile);

fprintf("Wrote: %s\n", outFile);