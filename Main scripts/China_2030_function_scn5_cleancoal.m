% Script file for the China electricity section modellling. This modelling
% takes into account the existing generation capacity (capacity as of 
% 2012 that would remain in operation in 2030.
%
% The script reads in the following input data:
% 1. Demand and RLDC for different PV and wind penetration from the MAT
%    file named "LDC_diffPV_wind_2030_China2.mat"
% 2. Random fuel and carbon prices and RLDC: 
%    "sample_fuel_carbon_price_mean_24.mat" and
%    "random_RLDC_xPV_xwind_full.mat"

%------------------------- created 09/07/13 -------------------------------
function China_2030_function_scn5_cleancoal(SD_IC)

%carbon_price_all = [0];
mean_carbon = 29;
NOX_price = 0;
SO2_price = 0;
PM_price = 0;

% NOX emission factors of each technology (g/MWh)
NOX_factor_coal = 478; 
NOX_factor_CCGT = 439; 
NOX_factor_nuke = 0;
NOX_factor_IGCC = 182;
NOX_factor_coal_exist = 478;
NOX_factor_CCGT_exist = 439; 
NOX_factor_nuke_exist = 0;

NOX_factor_with_exist = [NOX_factor_coal NOX_factor_CCGT NOX_factor_nuke ...
    NOX_factor_IGCC NOX_factor_coal_exist NOX_factor_CCGT_exist ...
    NOX_factor_nuke_exist];

% SO2 emission factors of each technology (g/MWh)
SO2_factor_coal = 649; 
SO2_factor_CCGT = 49; 
SO2_factor_nuke = 0;
SO2_factor_IGCC = 8;
SO2_factor_coal_exist = 649;
SO2_factor_CCGT_exist = 49; 
SO2_factor_nuke_exist = 0;

SO2_factor_with_exist = [SO2_factor_coal SO2_factor_CCGT SO2_factor_nuke ...
    SO2_factor_IGCC SO2_factor_coal_exist SO2_factor_CCGT_exist ...
    SO2_factor_nuke_exist];

% PM (2.5) emission factors of each technology (g/MWh)
PM_factor_coal = 53; 
PM_factor_CCGT = 0; 
PM_factor_nuke = 0;
PM_factor_IGCC = 22.7;
PM_factor_coal_exist = 53;
PM_factor_CCGT_exist = 0; 
PM_factor_nuke_exist = 0;

PM_factor_with_exist = [PM_factor_coal PM_factor_CCGT PM_factor_nuke ...
    PM_factor_IGCC PM_factor_coal_exist PM_factor_CCGT_exist ...
    PM_factor_nuke_exist];

% NOX, SOX and PM control costs for each fuel ($/MWh). This cost will be 
% included in the variable costs
extern_env_cost_coal = 4.67;
extern_env_cost_gas = 0;

% PV_wind_pen_all = [5 5; 5 10; 10 20; 20 30; 30 40];
PV_wind_pen_all = [30 40];

% for ii = 1:length(carbon_price_all)
    for jj = 1:size(PV_wind_pen_all,1)
%        mean_carbon = carbon_price_all(ii);
        PV_wind_pen = PV_wind_pen_all(jj,:);

%% ***********************************************************************%
% Specify input variables
%*************************************************************************
% Percentage uncertainty of carbon price and demand
SD_carbon_percent = 0.5; % based on the assumption
SD_demand_percent = 0.04; % SD of demand

% % User input for mean carbon price and wind and PV penetration level
% mean_carbon = input('Input mean carbon price ($/ton of CO2) = '); 
% PV_wind_pen = input('Input PV and wind penetration scenario = ');

SD_carbon = mean_carbon * SD_carbon_percent; 
SD_IC = 2.576; % number of SD to achieve 99% reliability criteria

% Input NOX, SO2 and PM costs

%% **********************************************************************
% Load random fuel and carbon prices, capital costs and RLDC for each PV
% and wind penetration scenario
%************************************************************************
% Specify path for data and sample files
LDC_data_path = ['D:\My Documents\Aust-China work\MATLAB\data\'];
sample_path = ['D:\My Documents\Aust-China work\MATLAB\Sample_file\'];

% load random fuel and carbon prices and capital costs
load([sample_path,['sample_fuel_carbon_price_mean_', ...
    num2str(mean_carbon),'_gasprice_9_coalprice_4.5.mat']])
load([sample_path,['sample_cap_cost_China.mat']])

% Load random RLDC 
vars_random_RLDC = {['random_RLDC_lesshydro_', num2str(PV_wind_pen(1)),...
    'PV_', num2str(PV_wind_pen(2)),'wind']};
load([sample_path,['random_RLDC_',num2str(PV_wind_pen(1)),'PV_', ...
    num2str(PV_wind_pen(2)),'wind_full.mat']],vars_random_RLDC{:});
eval(['random_RLDC', '= random_RLDC_lesshydro_' num2str(PV_wind_pen(1)) 'PV_' ...
    num2str(PV_wind_pen(2)) 'wind' ';']) % assign variable RLDC 
clear (['random_RLDC_lesshydro_',num2str(PV_wind_pen(1)),'PV_', ...
    num2str(PV_wind_pen(2)),'wind']);

%% ***********************************************************************
% Load LDC and RLDC for different wind and PV penetrations
% ************************************************************************
% Specify variables in the 'LDC_diffPV_wind_2030_China_2.mat' to be loaded
vars = {'Demand_hourly_2030','Demand_2030_less_existwind', ...
    'Demand_residual_diffscn','Demand_residual_15SC_chrono_diffscn', ...
    'Hydro_hourly_2030', ...
    'Energy_project_2030','Energy_PV_max_diffscn', ...
    'Energy_wind_max_diffscn','Energy_hydro_project_2030', ...
    'Energy_actual_PV_diffscn','Energy_actual_wind_diffscn', ...
    'Energy_actual_PVwind_diffscn', ...
    'Energy_actual_PVwind_lesshydro_diffscn', ... 
    'Energy_actual_newPV_diffscn','Energy_actual_existPV_diffscn', ...
    'Energy_actual_newwind_diffscn','Energy_actual_existwind_diffscn', ...
    'Energy_spill_diffscn_below15SC_lesshydro', ...
    'Energy_spill_diffscn_below0_lesshydro', ...
    'LDC_15pc_2030','LDC_85pc_2030','LDC_2030', ...
    'RLDC_2030_less_existwindPV', ...
    'IC_PV_existing','IC_PV_total_diffpen','IC_wind_total_diffpen', ...
    'IC_Wind_existing','IC_hydro_2030', ... 
    'PV_wind_gen_diffscn','PV_wind_gen_15SC_diffscn', ...
    'RLDC_diffscn','RLDC_diffscn_15SC', ...
    'RLDC_diffscn_15SC_lesshydro','RLDC_diffscn_lesshydro', ... 
    'PV_CF_desired','Wind_CF_desired','CF_hydro_actual_2012',...
    'percent_period_RLDC_below0_lesshydro_diffscn', ...
    'percent_period_RLDC_below15SC_lesshydro_diffscn', ...
    'x_axis_demand','x_axis_demand_percent'};

load([LDC_data_path,'LDC_diffPV_wind_2030_China_2.mat'], vars{:})

% Assign variables (for each PV and wind penetration scenario)
if PV_wind_pen == [5 5]
    RLDC_15SC = RLDC_diffscn_15SC(1,:);
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(1,:);
    Energy_actual_PVwind = Energy_actual_PVwind_diffscn(1);
    Energy_actual_PV = Energy_actual_PV_diffscn(1);
    Energy_actual_wind = Energy_actual_wind_diffscn(1);
    Energy_actual_existPV = Energy_actual_existPV_diffscn(1);
    Energy_actual_newPV = Energy_actual_newPV_diffscn(1);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(1);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(1);
    Energy_spill_PVwind = Energy_spill_diffscn_below15SC_lesshydro(1);
    IC_PV = IC_PV_total_diffpen(1);
    IC_Wind = IC_wind_total_diffpen(1); 
elseif PV_wind_pen == [5 10]
    RLDC_15SC = RLDC_diffscn_15SC(2,:);
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(2,:);
    Energy_actual_PVwind = Energy_actual_PVwind_diffscn(2);
    Energy_actual_PV = Energy_actual_PV_diffscn(2);
    Energy_actual_wind = Energy_actual_wind_diffscn(2);
    Energy_actual_existPV = Energy_actual_existPV_diffscn(2);
    Energy_actual_newPV = Energy_actual_newPV_diffscn(2);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(2);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(2);
    Energy_spill_PVwind = Energy_spill_diffscn_below15SC_lesshydro(2);
    IC_PV = IC_PV_total_diffpen(1);
    IC_Wind = IC_wind_total_diffpen(2);
elseif PV_wind_pen == [10 20]
    RLDC_15SC = RLDC_diffscn_15SC(3,:);
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(3,:);
    Energy_actual_PVwind = Energy_actual_PVwind_diffscn(3);
    Energy_actual_PV = Energy_actual_PV_diffscn(3);
    Energy_actual_wind = Energy_actual_wind_diffscn(3);
    Energy_actual_existPV = Energy_actual_existPV_diffscn(3);
    Energy_actual_newPV = Energy_actual_newPV_diffscn(3);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(3);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(3);
    Energy_spill_PVwind = Energy_spill_diffscn_below15SC_lesshydro(3);
    IC_PV = IC_PV_total_diffpen(2);
    IC_Wind = IC_wind_total_diffpen(3);
elseif PV_wind_pen == [20 30]
    RLDC_15SC = RLDC_diffscn_15SC(4,:);
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(4,:);
    Energy_actual_PVwind = Energy_actual_PVwind_diffscn(4);
    Energy_actual_PV = Energy_actual_PV_diffscn(4);
    Energy_actual_wind = Energy_actual_wind_diffscn(4);
    Energy_actual_existPV = Energy_actual_existPV_diffscn(4);
    Energy_actual_newPV = Energy_actual_newPV_diffscn(4);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(4);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(4);
    Energy_spill_PVwind = Energy_spill_diffscn_below15SC_lesshydro(4);
    IC_PV = IC_PV_total_diffpen(3);
    IC_Wind = IC_wind_total_diffpen(4);
elseif PV_wind_pen == [30 40]
    RLDC_15SC = RLDC_diffscn_15SC(5,:);
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(5,:);
    Energy_actual_PVwind = Energy_actual_PVwind_diffscn(5);
    Energy_actual_PV = Energy_actual_PV_diffscn(5);
    Energy_actual_wind = Energy_actual_wind_diffscn(5);
    Energy_actual_existPV = Energy_actual_existPV_diffscn(5);
    Energy_actual_newPV = Energy_actual_newPV_diffscn(5);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(5);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(5);
    Energy_spill_PVwind = Energy_spill_diffscn_below15SC_lesshydro(5);
    IC_PV = IC_PV_total_diffpen(4);
    IC_Wind = IC_wind_total_diffpen(5);
end
        
Energy_actual_wind_all = Energy_actual_wind; %(TWh)
Energy_actual_PV_all = Energy_actual_PV;
demand = RLDC_15SC_lesshydro;
demand_1 = demand;
hour = length(demand);

%% ***********************************************************************
% Define number of generation portfolio and increment of the share of each
% technology in the portfolio and LDC intervals
%*************************************************************************
% Four fossil-fuel generating technologies: Coal, CCGT, Nuclear and IGCC
gen_data = [];
for a = 0:10:100
    for b = 0:10:100
        for c = 0:10:100
            for d = 0:10:100
                total_mix = a+b+c+d;
                if total_mix == 100
                    gen_data = [gen_data; a b c d];
                end
            end
        end
    end
end

% Combined energy of PV, wind and hydro
Energy_PV_Wind_Hydro = Energy_actual_PV + Energy_actual_wind + ...
    Energy_hydro_project_2030; % (TWh)

% Specify number of random sample (number of Monte Carlo runs) and number
% of hours in each interval of the LDC (RLDCs)
interval = 1; % no. of interval in each bins in the average LDC and RLDCs
LDC_sample = 10000;
% LDC_sample = 5;
year_interval = 8760/interval; % number of LDC interval in a year
Ntech = 4; % no. of fossil-fuel technology considered

% Define minimum generation of each conventional technology
min_Coal = 0; min_CCGT = 0; min_Nuke = 0; min_IGCC = 0;
min_Coal_exist = 0; min_CCGT_exist = 0; min_Nuke_exist = 0;
min_Hydro = 0; 

min_gen = [min_Coal; min_CCGT; min_Nuke; min_IGCC];
min_gen_with_exist = [min_Coal; min_CCGT; min_Nuke; min_IGCC; ...
    min_Coal_exist; min_CCGT_exist; min_Nuke_exist;];

% Total energy for each sample of RLDC
Energy_random_RLDC = (sum(random_RLDC,2) * interval)/1000; % in TWh

%% **********************************************************************
% Existing generation fleet in 2030
%*************************************************************************
% Capacity that still exists in 2030 (GW)
IC_coal_exist = 676; IC_CCGT_exist = 37; IC_Nuke_exist = 13;
IC_IGCC_exist = 0; IC_Hydro_exist = 250; 

IC_exist = [IC_coal_exist IC_CCGT_exist IC_Nuke_exist IC_IGCC_exist];
%IC_exist_fossil = [IC_coal_exist IC_CCGT_exist IC_OCGT_exist];

IC_exist_port = repmat(IC_exist,length(gen_data),1);
%IC_exist_fossil_port = repmat(IC_exist_fossil,length(gen_data),1);

% New capacity of new PV, wind and Hydro 
IC_PV_new = IC_PV - IC_PV_existing;
IC_Wind_new = IC_Wind - IC_Wind_existing;
IC_Hydro_new = IC_hydro_2030 - IC_Hydro_exist;

%% ************************************************************************
% Installed generation capacity
%**************************************************************************
% Calculate the installed capacity of conventional generation technologies
% (CCGT and OCGT) to satisfy demand x% of the time.
mean_peak_demand = max(RLDC_15SC_lesshydro); %specify the peak demand (MW)
SD_demand = SD_demand_percent * mean_peak_demand; %specify the SD of peak demand (MW)

% Criteria for determining the installed capacity based on how many
% SD from the mean peak demand (from normal dist. table. SD_IC is the user 
% input and can be selected from the normal dist. table
%       90% reliability = 1.645 SD : demand is met 90% of the time
%       95% reliability = 1.96 SD
%       99% reliability = 2.576 SD
%       99.99% relia    = 3 SD

Z_SD = SD_IC; % criteria of determining installed capacity 

% Installed capacity is inclusive of all fossil-fuel technologies: coal, 
% gas, nuclear and IGCC (but without hydro)
IC_adjust = SD_demand * Z_SD;  
Installed_cap = mean_peak_demand + IC_adjust; % total 

% IC with hydro, PV and Wind
Installed_cap_withhydro = Installed_cap + IC_hydro_2030;
Installed_cap_total = Installed_cap_withhydro + IC_PV + IC_Wind;

% IC of each fossil fuel technology in each portfolio
Installed_cap_port = (gen_data/100) * Installed_cap;

% Real build capacity in 2030 (capacity that incurs capital costs: Coal,
% CCGT, Nuclear, IGCC)
IC_realbuild_port = Installed_cap_port - IC_exist_port;
IC_realbuild_port(IC_realbuild_port < 0)=0;

% Portfolios which exclude those with nuclear capacity is greater 
% than 198 GW (total new and existing)
port_exclude = find(IC_realbuild_port(:,3) > (0.5*Installed_cap_total)); 
IC_realbuild_port(port_exclude,:) = [];
IC_exist_port(port_exclude,:) = [];
Installed_cap_port(port_exclude,:) = [];

% Installed capacity of each technology categorized into new and existing
IC_port_with_exist = [IC_realbuild_port ...
    Installed_cap_port(:,1:3) - IC_realbuild_port(:,1:3)];
IC_port_with_exist_withhydro = [IC_port_with_exist ...
    repmat(IC_hydro_2030,size(IC_port_with_exist,1),1)];
IC_port_with_exist_withPVWind_hydro = [IC_port_with_exist ...
    repmat(IC_PV,size(IC_port_with_exist,1),1) ...
    repmat(IC_Wind,size(IC_port_with_exist,1),1) ...
    repmat(IC_hydro_2030,size(IC_port_with_exist,1),1)];

% Share of technology categorised into new and existing
IC_percent_port_withexist = ...
    (IC_port_with_exist/Installed_cap)*100;
IC_percent_port_withexist_total = ...
    (IC_port_with_exist_withPVWind_hydro/Installed_cap_total)*100;

% IC with every technology (combining existing and new)
Installed_cap_port_total = [Installed_cap_port ...
    IC_port_with_exist_withPVWind_hydro(:,8:10)];
IC_percent_port_total = ...
    (Installed_cap_port_total/Installed_cap_total)*100;

%% ************************************************************************
% Technical parameters of each generation technology
%*************************************************************************
% Overnight capital cost ($/MW)
mean_cap_cost_coal = 511000;
mean_cap_cost_CCGT = 391000; 
mean_cap_cost_nuke = 3200000;
mean_cap_cost_IGCC = 1206000;
mean_cap_cost_PV = 830000;
mean_cap_cost_wind = 530000;
mean_cap_cost_hydro = 982000;

% Fixed O&M cost of each technology ($/MW/year)
FOM_coal = 29361; FOM_CCGT = 41036; FOM_nuke = 153630; FOM_IGCC = 78365; 
FOM_coal_exist = 30448; FOM_CCGT_exist = 43013; FOM_nuke_exist = 153630;
FOM_PV = 14916; FOM_wind = 38115; FOM_hydro = 49102;

FOM_with_exist = [FOM_coal FOM_CCGT FOM_nuke FOM_IGCC FOM_coal_exist ...
    FOM_CCGT_exist FOM_nuke_exist];

% Plant life of each technology (years)
plant_life_coal = 40; plant_life_CCGT = 30; plant_life_nuke = 40;
plant_life_IGCC = 20; plant_life_PV = 20; plant_life_wind = 30;
plant_life_hydro = 40;

% Variable O&M cost ($/MWh)
VOM_coal = 2.066; VOM_CCGT = 1.477; VOM_nuke = 2.216; VOM_IGCC = 2.954;
VOM_coal_exist = 2.175; VOM_CCGT_exist = 1.477; VOM_nuke_exist = 2.216; 
VOM_PV = 7.386; VOM_wind = 7.386; VOM_hydro = 8.863; 

VOM_fossil = [VOM_coal VOM_CCGT VOM_nuke VOM_IGCC];
VOM_fossil_with_exist = [VOM_fossil VOM_coal_exist VOM_CCGT_exist ...
    VOM_nuke_exist];

% Heat rate of each technology (GJ/MWh)
HR_coal = 9; HR_CCGT = 6.00; HR_nuke = 9.73; HR_IGCC = 7.83;
HR_coal_exist = 9.05; HR_CCGT_exist = 6.43; HR_nuke_exist = 9.83; 

Heat_rate = [HR_coal HR_CCGT HR_nuke HR_IGCC HR_coal_exist ...
    HR_CCGT_exist HR_nuke_exist];

% CO2 emission factors of each technology (tCO2/MWh)
CO2_factor_coal = 0.918; 
CO2_factor_CCGT = 0.382; 
CO2_factor_nuke = 0;
CO2_factor_IGCC = 0.896;
CO2_factor_coal_exist = 0.982;
CO2_factor_CCGT_exist = 0.411; 
CO2_factor_nuke_exist = 0;
CO2_factor_Hydro = 0;

CO2_factor_with_exist = [CO2_factor_coal CO2_factor_CCGT CO2_factor_nuke ...
    CO2_factor_IGCC CO2_factor_coal_exist CO2_factor_CCGT_exist ...
    CO2_factor_nuke_exist];

%% ************************************************************************
% Calculate annualised fixed cost of each technology
%*************************************************************************
n = LDC_sample;
WACC = 0.1; % Weighted Average Cost of Capital (WACC)

% Calculate the Capital Cost Recovery Factor (CRF)
CRF_coal = (WACC*((1+WACC)^plant_life_coal))/...
    (((1+WACC)^plant_life_coal)-1);
CRF_CCGT = (WACC*((1+WACC)^plant_life_CCGT))/...
    (((1+WACC)^plant_life_CCGT)-1);
CRF_nuke = (WACC*((1+WACC)^plant_life_nuke))/...
    (((1+WACC)^plant_life_nuke)-1);
CRF_IGCC = (WACC*((1+WACC)^plant_life_IGCC))/...
    (((1+WACC)^plant_life_IGCC)-1);
CRF_PV = (WACC*((1+WACC)^plant_life_PV))/...
    (((1+WACC)^plant_life_PV)-1);
CRF_wind = (WACC*((1+WACC)^plant_life_wind))/...
    (((1+WACC)^plant_life_wind)-1);
CRF_hydro = (WACC*((1+WACC)^plant_life_hydro))/...
    (((1+WACC)^plant_life_hydro)-1);

% Calculate annualized captial cost of each tech for the reference cost
annualized_cap_cost_coal_ref = mean_cap_cost_coal * CRF_coal;
annualized_cap_cost_CCGT_ref = mean_cap_cost_CCGT * CRF_CCGT;
annualized_cap_cost_nuke_ref = mean_cap_cost_nuke * CRF_nuke;
annualized_cap_cost_IGCC_ref = mean_cap_cost_IGCC * CRF_IGCC;
annualized_cap_cost_PV_ref = mean_cap_cost_PV * CRF_PV;
annualized_cap_cost_wind_ref = mean_cap_cost_wind * CRF_wind;
annualized_cap_cost_hydro_ref = mean_cap_cost_hydro * CRF_hydro;

annualized_cap_cost_ref = [annualized_cap_cost_coal_ref ...
    annualized_cap_cost_CCGT_ref annualized_cap_cost_nuke_ref ...
    annualized_cap_cost_PV_ref annualized_cap_cost_wind_ref ...
    annualized_cap_cost_hydro_ref];

% Call random overnight capital cost ($/MW) of each technology which have
% been generated and stored in a MAT file name random_capital_cost.mat 
% with the matrix name random_cap_cost_all
cap_cost_coal_all = random_cap_cost_all(:,1);
cap_cost_CCGT_all = random_cap_cost_all(:,2);
cap_cost_nuke_all = random_cap_cost_all(:,3);
cap_cost_IGCC_all = random_cap_cost_all(:,4);
cap_cost_PV_all = random_cap_cost_all(:,5);
cap_cost_wind_all = random_cap_cost_all(:,6);
cap_cost_hydro_all = random_cap_cost_all(:,7);

% Determine the SD of the random capital cost
SD_cap_cost_coal = std(cap_cost_coal_all);
SD_cap_cost_CCGT = std(cap_cost_CCGT_all);
SD_cap_cost_nuke = std(cap_cost_nuke_all);
SD_cap_cost_IGCC = std(cap_cost_IGCC_all);
SD_cap_cost_PV = std(cap_cost_PV_all);
SD_cap_cost_wind = std(cap_cost_wind_all);
SD_cap_cost_hydro = std(cap_cost_hydro_all);

% Annualized capital cost for every random capital cost ($/MW)
annualized_cap_cost_coal = (cap_cost_coal_all * CRF_coal);
annualized_cap_cost_CCGT = (cap_cost_CCGT_all * CRF_CCGT);
annualized_cap_cost_nuke = (cap_cost_nuke_all * CRF_nuke);
annualized_cap_cost_IGCC = (cap_cost_IGCC_all * CRF_IGCC);
annualized_cap_cost_PV = (cap_cost_PV_all * CRF_PV);
annualized_cap_cost_wind = (cap_cost_wind_all * CRF_wind);
annualized_cap_cost_hydro = (cap_cost_hydro_all * CRF_hydro);

random_annual_cap_cost_fossil = [annualized_cap_cost_coal ...
    annualized_cap_cost_CCGT annualized_cap_cost_nuke ...
    annualized_cap_cost_IGCC];
random_annual_cap_cost_all = [random_annual_cap_cost_fossil ...
    annualized_cap_cost_PV annualized_cap_cost_wind ...
    annualized_cap_cost_hydro];

%% ***********************************************************************
% Calculate generation cost based on the expected coal,gas and carbon 
% prices under the expected Load Duration Curve
%*************************************************************************
% Mean and SD of fuel price ($/GJ)
mean_coal = 4.52;
mean_gas = 8.938;
mean_nuke = 0.7;

% Calculate fuel cost of each tech for the reference case
fuel_cost_coal_ref = mean_coal * HR_coal; % ($/MWh)
fuel_cost_CCGT_ref = mean_gas * HR_CCGT;
fuel_cost_nuke_ref = mean_nuke * HR_nuke;
fuel_cost_IGCC_ref = mean_coal * HR_IGCC;
fuel_cost_coal_exist_ref = mean_coal * HR_coal_exist;
fuel_cost_CCGT_exist_ref = mean_gas * HR_CCGT_exist;
fuel_cost_nuke_exist_ref = mean_nuke * HR_nuke_exist;

fuel_cost_with_exist_ref = [fuel_cost_coal_ref fuel_cost_CCGT_ref ...
     fuel_cost_nuke_ref fuel_cost_IGCC_ref ...
     fuel_cost_coal_exist_ref fuel_cost_CCGT_ref fuel_cost_nuke_ref];

% Calculate carbon cost of each tech for the reference case
carbon_cost_coal_ref = mean_carbon * CO2_factor_coal;
carbon_cost_CCGT_ref = mean_carbon * CO2_factor_CCGT;
carbon_cost_nuke_ref = mean_carbon * CO2_factor_nuke;
carbon_cost_IGCC_ref = mean_carbon * CO2_factor_IGCC;
carbon_cost_coal_exist_ref = mean_carbon * CO2_factor_coal_exist;
carbon_cost_CCGT_exist_ref = mean_carbon * CO2_factor_CCGT_exist;
carbon_cost_nuke_exist_ref = mean_carbon * CO2_factor_nuke_exist;

carbon_cost_with_exist_ref = [carbon_cost_coal_ref ...
    carbon_cost_CCGT_ref carbon_cost_nuke_ref carbon_cost_IGCC_ref...
    carbon_cost_coal_exist_ref carbon_cost_CCGT_exist_ref ...
    carbon_cost_nuke_exist_ref];

% Calculate environmental damangd cost of each tech for the reference case
extern_env_cost_coal_ref = extern_env_cost_coal;
extern_env_cost_CCGT_ref = extern_env_cost_gas;
extern_env_cost_nuke_ref = 0;
extern_env_cost_IGCC_ref = extern_env_cost_coal;
extern_env_cost_coal_exist_ref = extern_env_cost_coal;
extern_env_cost_CCGT_exist_ref = extern_env_cost_gas;
extern_env_cost_nuke_exist_ref = 0;

extern_env_cost_with_exist_ref = [extern_env_cost_coal_ref ...
    extern_env_cost_CCGT_ref extern_env_cost_nuke_ref ...
    extern_env_cost_IGCC_ref extern_env_cost_coal_exist_ref ...
    extern_env_cost_CCGT_exist_ref extern_env_cost_nuke_exist_ref];

% Air Pollutant taxes
NOX_cost_coal_ref = NOX_price * (NOX_factor_coal/10^6);
NOX_cost_CCGT_ref = NOX_price * (NOX_factor_CCGT/10^6);
NOX_cost_nuke_ref = NOX_price * (NOX_factor_nuke/10^6);
NOX_cost_IGCC_ref = NOX_price * (NOX_factor_IGCC/10^6);
NOX_cost_coal_exist_ref = NOX_price * (NOX_factor_coal_exist/10^6);
NOX_cost_CCGT_exist_ref = NOX_price * (NOX_factor_CCGT_exist/10^6);
NOX_cost_nuke_exist_ref = NOX_price * (NOX_factor_nuke_exist/10^6);

NOX_cost_with_exist_ref = [NOX_cost_coal_ref ...
    NOX_cost_CCGT_ref NOX_cost_nuke_ref NOX_cost_IGCC_ref...
    NOX_cost_coal_exist_ref NOX_cost_CCGT_exist_ref ...
    NOX_cost_nuke_exist_ref];

SO2_cost_coal_ref = SO2_price * (SO2_factor_coal/10^6);
SO2_cost_CCGT_ref = SO2_price * (SO2_factor_CCGT/10^6);
SO2_cost_nuke_ref = SO2_price * (SO2_factor_nuke/10^6);
SO2_cost_IGCC_ref = SO2_price * (SO2_factor_IGCC/10^6);
SO2_cost_coal_exist_ref = SO2_price * (SO2_factor_coal_exist/10^6);
SO2_cost_CCGT_exist_ref = SO2_price * (SO2_factor_CCGT_exist/10^6);
SO2_cost_nuke_exist_ref = SO2_price * (SO2_factor_nuke_exist/10^6);

SO2_cost_with_exist_ref = [SO2_cost_coal_ref ...
    SO2_cost_CCGT_ref SO2_cost_nuke_ref SO2_cost_IGCC_ref...
    SO2_cost_coal_exist_ref SO2_cost_CCGT_exist_ref ...
    SO2_cost_nuke_exist_ref];

PM_cost_coal_ref = PM_price * (PM_factor_coal/10^6);
PM_cost_CCGT_ref = PM_price * (PM_factor_CCGT/10^6);
PM_cost_nuke_ref = PM_price * (PM_factor_nuke/10^6);
PM_cost_IGCC_ref = PM_price * (PM_factor_IGCC/10^6);
PM_cost_coal_exist_ref = PM_price * (PM_factor_coal_exist/10^6);
PM_cost_CCGT_exist_ref = PM_price * (PM_factor_CCGT_exist/10^6);
PM_cost_nuke_exist_ref = PM_price * (PM_factor_nuke_exist/10^6);

PM_cost_with_exist_ref = [PM_cost_coal_ref ...
    PM_cost_CCGT_ref PM_cost_nuke_ref PM_cost_IGCC_ref...
    PM_cost_coal_exist_ref PM_cost_CCGT_exist_ref ...
    PM_cost_nuke_exist_ref];

% Calculate Variable cost for the reference case ($/MWh) without carbon
% price
VC_pu_fossil_with_exist_ref = VOM_fossil_with_exist + ...
    fuel_cost_with_exist_ref + extern_env_cost_with_exist_ref +...
    NOX_cost_with_exist_ref + SO2_cost_with_exist_ref + ...
    PM_cost_with_exist_ref;

%% *********************************************************************** 
% Use the correlated random fuel and carbon prices which have been
% generated and stored in a MAT file name sample_fuel_carbon_price_meanx.mat
% with the matrix name sample_fuel_carbon_price
%*************************************************************************
% Specify fuel and carbon price from sample_fuel_carbon_price matrix
coal_price = sample_fuel_carbon_price(:,1); % ($/GJ)
gas_price = sample_fuel_carbon_price(:,2);
nuke_price = sample_fuel_carbon_price(:,3);
carbon_price = sample_fuel_carbon_price(:,4);

% Calculate fuel cost of each technolgy for every sample ($/MWh)
fuel_cost_coal = coal_price * HR_coal;
fuel_cost_CCGT = gas_price * HR_CCGT;
fuel_cost_nuke = nuke_price * HR_nuke;
fuel_cost_IGCC = coal_price * HR_IGCC;
fuel_cost_coal_exist = coal_price * HR_coal_exist;
fuel_cost_CCGT_exist = gas_price * HR_CCGT_exist;
fuel_cost_nuke_exist = nuke_price * HR_nuke_exist;
fuel_cost_Hydro = zeros(LDC_sample,1);

fuel_cost_with_exist = [fuel_cost_coal fuel_cost_CCGT ...
    fuel_cost_nuke fuel_cost_IGCC ...
    fuel_cost_coal_exist fuel_cost_CCGT_exist fuel_cost_nuke_exist];

% Calculate carbon cost of each technology for every sample ($/MWh)
carbon_cost_coal = carbon_price * CO2_factor_coal;
carbon_cost_CCGT = carbon_price * CO2_factor_CCGT;
carbon_cost_nuke = carbon_price * CO2_factor_nuke;
carbon_cost_IGCC = carbon_price * CO2_factor_IGCC;
carbon_cost_coal_exist = carbon_price * CO2_factor_coal_exist;
carbon_cost_CCGT_exist = carbon_price * CO2_factor_CCGT_exist;
carbon_cost_nuke_exist = carbon_price * CO2_factor_nuke_exist;
carbon_cost_Hydro = zeros(LDC_sample,1);

carbon_cost_with_exist = [carbon_cost_coal_exist carbon_cost_CCGT ...
    carbon_cost_nuke carbon_cost_IGCC carbon_cost_coal_exist ...
    carbon_cost_CCGT_exist carbon_cost_nuke_exist ];

% calculate variable operating cost per unit ($/MWh) of each tech. for 
% every sample ($/MWh)
VOM_all_fossil_with_exist = ones(LDC_sample,1) * VOM_fossil_with_exist;
VC_pu_with_exist = VOM_all_fossil_with_exist + fuel_cost_with_exist + ...
    carbon_cost_with_exist + repmat(extern_env_cost_with_exist_ref,n,1) + ...
    repmat(NOX_cost_with_exist_ref,n,1) + ...
    repmat(SO2_cost_with_exist_ref,n,1) + ...
    repmat(PM_cost_with_exist_ref,n,1);

% Variable operating cost of each technology for every set of samples
% ($/MWh)
VC_pu_coal = VC_pu_with_exist(:,1);
VC_pu_CCGT = VC_pu_with_exist(:,2); 
VC_pu_nuke = VC_pu_with_exist(:,3); 
VC_pu_IGCC = VC_pu_with_exist(:,4); 
VC_pu_coal_exist = VC_pu_with_exist(:,5);
VC_pu_CCGT_exist = VC_pu_with_exist(:,6);
VC_pu_nuke_exist = VC_pu_with_exist(:,7); 

%% *********************************************************************** 
% Probability that demand is met when using the adjusted LDC 
%**************************************************************************
prob_demand_met = 1-(sum(random_RLDC(:,1) > Installed_cap)/LDC_sample);
expected_peak_demand = (sum(random_RLDC(:,1)))/LDC_sample;
SD_random_peak = std(random_RLDC(:,1));
SD_random_peak_percent = SD_random_peak/expected_peak_demand; 

%% ***********************************************************************
% Fixed and variable cost of PV, wind and Hydro generation (not included in
% the dispatch)
%*************************************************************************
% Total variable costs of PV, wind and hydro (dollars)
total_VC_PV_new = (Energy_actual_newPV * 10^6 * VOM_PV); % $
total_VC_PV_exist = (Energy_actual_existPV * 10^6 * VOM_PV);
total_VC_PV_all = total_VC_PV_new + total_VC_PV_exist;
total_VC_PV_all_everyport = repmat(total_VC_PV_all, ...
    length(IC_port_with_exist),1);

total_VC_wind_new = (Energy_actual_newwind * 10^6 * VOM_wind); % $
total_VC_wind_exist = (Energy_actual_existwind * 10^6 * VOM_wind);
total_VC_wind_all = total_VC_wind_new + total_VC_wind_exist;
total_VC_wind_all_everyport = repmat(total_VC_wind_all, ...
    length(IC_port_with_exist),1);

total_VC_hydro = (Energy_hydro_project_2030 * 10^6 * VOM_hydro); % $
total_VC_hydro_everyport = repmat(total_VC_hydro, ...
    length(IC_port_with_exist),1);

% Combined total variable costs ($) of PV, wind and hydro
total_VC_PV_Wind_Hydro = total_VC_PV_all + total_VC_wind_all + ...
    total_VC_hydro;
total_VC_PV_Wind_Hydro_everyport = repmat(total_VC_PV_Wind_Hydro, ...
    length(IC_port_with_exist),1);
total_VC_PV_Wind_Hydro_everysample = ...
    repmat(total_VC_PV_Wind_Hydro,LDC_sample,1);

% Total random FC for PV, wind and hydro ($) which calculate for every
% random annualized capital costs (This is the same for every generation
% portfolio)
total_FOM_PV_new = FOM_PV * IC_PV_new * 10^3; %($)
total_FOM_PV_existing = FOM_PV * IC_PV_existing * 10^3;
total_capcost_PV_new_everysample = ...
    repmat((IC_PV_new * 10^3),LDC_sample,1).* ...
    random_annual_cap_cost_all(:,5);
total_FC_PV_all_everysample = total_capcost_PV_new_everysample + ...
    total_FOM_PV_existing + total_FOM_PV_new;
    
total_FOM_wind_new = FOM_wind * IC_Wind_new * 10^3; %($)
total_FOM_wind_existing = FOM_wind * IC_Wind_existing * 10^3;
total_capcost_wind_new_everysample = ...
    repmat((IC_Wind_new * 10^3),LDC_sample,1).* ...
    random_annual_cap_cost_all(:,6);
total_FC_wind_all_everysample = total_capcost_wind_new_everysample + ...
    total_FOM_wind_existing + total_FOM_wind_new ;

total_FOM_hydro_new = FOM_hydro * IC_Hydro_new * 10^3; %($)
total_FOM_hydro_existing = FOM_wind * IC_Hydro_exist * 10^3;
total_capcost_hydro_new_everysample = ...
    repmat((IC_Hydro_new * 10^3),LDC_sample,1).* ...
    random_annual_cap_cost_all(:,7);
total_FC_hydro_all_everysample = total_capcost_hydro_new_everysample + ...
    total_FOM_hydro_existing + total_FOM_hydro_new;

total_FC_PV_Wind_Hydro_everysample = total_FC_PV_all_everysample + ...
    total_FC_wind_all_everysample + total_FC_hydro_all_everysample;

%% ************************************************************************
% Determine merit order dispatch of each technology based on their VC
%*************************************************************************
tech = {'coal','CCGT','nuke','IGCC','existing coal','existing CCGT', ...
    'existing nuke'};    
cost_ENS = 1000; % unit cost of energy not served ($/MWh)

tic
for r = 1:length(IC_exist_port) % For every generation portfolio 
IC_with_exist = IC_port_with_exist(r,:);
    IC_realbuild = IC_realbuild_port(r,:);
% for r = 257:258 
%     IC_with_exist = IC_port_with_exist(258,:);
%     IC_realbuild = IC_realbuild_port(258,:);
    cap_tech = IC_with_exist - min_gen_with_exist'; %available capacity of each technology
    gen_output = zeros(length(cap_tech),length(demand)); %initialize matrix
    total_capcost_fossil = repmat((IC_realbuild * 10^3),LDC_sample,1).* ...
        random_annual_cap_cost_fossil;
    total_FOM_fossil = IC_port_with_exist(r,:) * 10^3 .* FOM_with_exist; % in $
    total_FC_fossil = [total_capcost_fossil zeros(LDC_sample,3)] + ...
        repmat(total_FOM_fossil,LDC_sample,1); % in $
    total_FC = [total_FC_fossil total_FC_PV_Wind_Hydro_everysample]; % for every sample

    for k = 1:size(random_RLDC,1) %For each set of random parameters (fuel, carbon prices and LDC)
        LDC = random_RLDC(k,:);
        LDC_1 = LDC;
        [VC_sort,idx] = sort(VC_pu_with_exist(k,:)); % Rank each technology based on VC
        cap_sort = cap_tech(idx); % capacity of each tech after sorted based on VC
        cap(k,:) = cap_sort;
        
        for i=1:length(cap_tech) % Determine output(GW) of each technology in each segment of the LDC
            gen_output(i,:) = min(LDC,cap_sort(i).*ones(1,hour));
            LDC = LDC - gen_output(i,:);
            [sortidx, ind] = sort(idx); % sort the index to arrange back to the original list
            gen_tech = gen_output(ind,:); % output (GW) in the order coal->CCGT->Nuke->IGCC
        end
        
        total_gen = sum(gen_tech,1); % in GW
        ENS_period = sum((total_gen*1.000000000009) < LDC_1); % Determine the number of interval that load is not met
        ENS_period_all(k,:) = ENS_period;
        ENS = (LDC_1(1:ENS_period) - Installed_cap)*interval; % amount of energy not served in a year (GWh) which is the same for every portfolio
        ENS_all(k,:) = sum(ENS); % in GWh
        ENS_cost = sum(ENS * cost_ENS * 10^3); % total cost of unserved energy ($)
        ENS_cost_all(k,:) = (ENS_cost/1000000); % ($Million)
        energy_tech_fossil = gen_tech * interval; % energy of each tech in each segment of each  LDC (GWh)
        energy_tech_year_fossil = sum(energy_tech_fossil,2); % energy of each technology in a year for a set of random parameters (GWh)
        energy_fossil(k,:) = energy_tech_year_fossil; % annual energy of each tech. for every sample of fuel and carbon price (GWh)
        VC_fossil = (energy_tech_year_fossil'.*VC_pu_with_exist(k,:))*10^3; % VC of each tech for each carbon price
        VC_tech_fossil(k,:) = VC_fossil; % in $
        CF = (energy_tech_year_fossil'./(IC_with_exist * length(LDC_1)* ...
            interval)) * 100;
        CF_all(k,:) = CF; % CF of each tech for every set of random parameters
        CF_all(isnan(CF_all))=0;
        CF_mean = mean(CF_all,1); % Mean CF for a generation portfolio
        LDC = LDC_1;
    end
    CF_mean_port(r,:) = CF_mean;
    VC_tech_withPV_Wind_Hydro = [VC_tech_fossil ...
          total_VC_PV_Wind_Hydro_everysample]; % VC ($) for every set of random parameters   
    VC_mean_port_fossil(r,:) = sum(VC_tech_fossil,1)/LDC_sample; % mean VC by each technology for every generation portfolio
    FC_mean_port_fossil(r,:) = sum(total_FC_fossil,1)/LDC_sample; % in $
%     VC_mean_port_fossil = sum(VC_tech_fossil,1)/LDC_sample; % mean VC by each technology for every generation portfolio
%     FC_mean_port_fossil = sum(total_FC_fossil,1)/LDC_sample; % in $
    totalcost_fossil = VC_tech_fossil + total_FC_fossil; % in $
    totalcost_withPV_Wind_Hydro = VC_tech_withPV_Wind_Hydro + total_FC; % in $
    totalcost_tech_fossil = totalcost_fossil/1000000; % total cost of each technology ($Million) 
    totalcost_tech_withPV_Wind_Hydro = totalcost_withPV_Wind_Hydro/1000000; % $ Million
    mean_energy_fossil(r,:) = sum(energy_fossil,1)/LDC_sample; % Mean energy for every generation portfolio (GWh)
%     mean_energy_fossil = sum(energy_fossil,1)/LDC_sample;    
    CO2_tech = ((repmat(CO2_factor_with_exist,[LDC_sample,1])).* ... 
        energy_fossil) * 10^3; % in ton of CO2
    total_CO2 = sum(CO2_tech,2); % total CO2 emission of generation portfolio (ton of CO2)
    NOX_tech = ((repmat(NOX_factor_with_exist,[LDC_sample,1])).* ...
        energy_fossil) * 10^3; % in gram of NOX 
    total_NOX = sum(NOX_tech,2);
    SO2_tech = ((repmat(SO2_factor_with_exist,[LDC_sample,1])).* ...
        energy_fossil) * 10^3; % in gram of SO2
    total_SO2 = sum(SO2_tech,2);
    PM_tech = ((repmat(PM_factor_with_exist,[LDC_sample,1])).* ...
        energy_fossil) * 10^3; % in ton of PM
    total_PM = sum(PM_tech,2);
    overallcost_mix_fossil(r,:) = sum(totalcost_tech_fossil,2) + ...
        ENS_cost_all; % overall generation cost ($Million) for every set of uncertain parameters (include cost of unserved energy)
    overallcost_mix_withPV_Wind_Hydro(r,:) = ...
        sum(totalcost_tech_withPV_Wind_Hydro,2) + ENS_cost_all; % in $Million
%     overallcost_mix_fossil = sum(totalcost_tech_fossil,2) + ...
%         ENS_cost_all; 
end
toc
overallcost_mix_billion = overallcost_mix_withPV_Wind_Hydro/1000;

% Total average variable cost and fixed cost by technology for every 
% generation portfolio ($ Million)
VC_mean_port = [VC_mean_port_fossil total_VC_PV_Wind_Hydro_everyport(1:r)];
avg_VC_port = sum(VC_mean_port,2)/1000000; %($ Million)

avg_FC_PV_Wind_Hydro = sum(total_FC_PV_Wind_Hydro_everysample)/LDC_sample; % average fixed cost is the same for every generation portfolio
FC_port = [FC_mean_port_fossil (avg_FC_PV_Wind_Hydro.*ones(r,1))]; % for every portfolio
avg_FC_port = sum(FC_port,2)/1000000; % ($ Million)

% Average total cost of every generation portfolio. 
%(***** This value is the same as cost_mean_mix *****)
avg_FC_VC_port = [avg_FC_port avg_VC_port];
avg_totalcost = sum(avg_FC_VC_port,2); 

% Mean energy (TWh) of each technology for each generation portfolio 
mean_energy_withPV_Wind_Hydro = [mean_energy_fossil/10^3 ...
    (Energy_PV_Wind_Hydro * ones(size(mean_energy_fossil,1),1))]; % TWh
mean_energy_total = sum(mean_energy_withPV_Wind_Hydro,2); % TWh

% Total energy demand (TWh) to be met by fossil-fuel gen and others in a year
total_energy_fossil = (sum(demand) * interval)/10^3; % TWh
total_energy = total_energy_fossil + Energy_PV_Wind_Hydro; % TWh

% Expected cost and SD of every generation portfolio ($Million)
cost_mean_mix = sum(overallcost_mix_withPV_Wind_Hydro,2)/LDC_sample; 
SD_mix = std(overallcost_mix_withPV_Wind_Hydro,0,2);

% Expected cost and of every generation portfolio expressed in $Billion
cost_mean_mix_billion = cost_mean_mix/1000;
SD_mix_billion = SD_mix/1000; 

% proportion of the cost of unserved energy to total cost
mean_ENS_cost = sum(ENS_cost_all)/length(random_RLDC); %Mean cost of energy not served which is the same for every generation portfolio
mean_ENS_cost_all = mean_ENS_cost*ones(r,1);
proportion_ENS_cost = ((mean_ENS_cost_all/1000000)./cost_mean_mix)*100; %Percentage of EN

%% ***********************************************************************
% Expected cost in per unit ($/MWh) of each generation portfolio
%*************************************************************************
energy_all_random_RLDC = Energy_random_RLDC + Energy_actual_PV + ...
    Energy_actual_wind + Energy_hydro_project_2030; % in TWh
overallcost_mix_pu = (overallcost_mix_withPV_Wind_Hydro*10^6)./ ...
    repmat((energy_all_random_RLDC * 10^6),1,length(IC_exist_port))';
% overallcost_mix_pu = (overallcost_mix_withPV_Wind_Hydro(258,:)*10^6)./ ...
%     repmat((energy_all_random_RLDC* 10^6),1,1)';
cost_mean_mix_pu = sum(overallcost_mix_pu,2)/LDC_sample; % ($/MWh)
SD_pu = std(overallcost_mix_pu,0,2);

%% ************************************************************************
% Calculate CO2, NOX, SO2 emission and PM of each generation portfolio
%*************************************************************************
for u = 1:length(mean_energy_fossil)
    CO2(u,:) = (CO2_factor_with_exist .* mean_energy_fossil(u,:)) * 10^3;
    NOX(u,:) = (NOX_factor_with_exist .* mean_energy_fossil(u,:)) * 10^3; % g
    SO2(u,:) = (SO2_factor_with_exist .* mean_energy_fossil(u,:)) * 10^3;% g
    PM(u,:) = (PM_factor_with_exist .* mean_energy_fossil(u,:)) * 10^3; % ton
end

totalCO2 = sum(CO2,2)/1000000; % total CO2 emissions (M.ton) for each generation mix
totalNOX = sum(NOX,2)/10^12; % in Mton 
totalSO2 = sum(SO2,2)/10^12; % in Mton
totalPM = sum(PM,2)/10^12; % in Mton

%% ***********************************************************************
% Statistical parameters
% ************************************************************************
%Generation cost at the 95th percentile
percentile_95 = prctile(overallcost_mix_billion,95,2);
[sort_percentile_95,port_sort_index_percentile] = ...
    sort(percentile_95,'descend');

% Skewness and kurtosis
skew_cost_mix = skewness(overallcost_mix_billion,0,2);
kurtosis_cost_mix = kurtosis(overallcost_mix_billion,0,2);

%% *********************************************************************** 
% Rearrange the results
% ************************************************************************
% sort the cost of each generation portfolios
[cost_mix_sort,port_cost_index] = sort(cost_mean_mix_billion,'descend');
[SD_mix_sort,port_SD_index] = sort(SD_mix_billion,'descend');
[cost_mix_pu_sort,port_cost_pu_index] = sort(cost_mean_mix_pu,'descend');
[SD_mix_pu_sort,port_SD_pu_index] = sort(SD_pu,'descend');

% rank portfolios in terms of cost and SD and display them in % 
port_cost_rank = IC_port_with_exist(port_cost_index,:);
port_SD_rank = gen_data(port_SD_index,:);
port_cost_pu_rank = IC_port_with_exist(port_cost_pu_index,:);
port_SD_pu_rank = gen_data(port_SD_pu_index,:);

% sort portfolio from high to low
port_cost_SD_CO2_NOX_SO2_PM_idx = [cost_mean_mix_billion(port_cost_index) ...
    SD_mix_billion(port_cost_index) totalCO2(port_cost_index) ...
    totalNOX(port_cost_index) totalSO2(port_cost_index) ...
    totalPM(port_cost_index) port_cost_index];
port_cost_SD_CO2_NOX_SO2_PM_pu_idx = [cost_mean_mix_pu(port_cost_pu_index) ...
    SD_pu(port_cost_pu_index) totalCO2(port_cost_pu_index) ...
    totalNOX(port_cost_pu_index) totalSO2(port_cost_pu_index) ...
    totalPM(port_cost_pu_index) port_cost_pu_index];

% matrix containing cost, SD, CO2 of the portfolio
cost_SD_CO2_NOX_SO2_PM_idx_port = [cost_mean_mix_billion SD_mix_billion ...
    totalCO2 totalNOX totalSO2 totalPM];
cost_SD_CO2_NOX_SO2_PM_idx_pu_port = [cost_mean_mix_pu SD_pu ...
    totalCO2 totalNOX totalSO2 totalPM];

%% ***********************************************************************
% Identify generation portfolios on the cost-risk efficient frontier
%*************************************************************************
% calling the function "prtp"  to identify portfolios along the efficient
cost_SD_mix_billion = [cost_mean_mix_billion SD_mix_billion ];
[cost_SD_port_EF_temp, idx_port_EF_temp] = prtp(cost_SD_mix_billion);  

% Put cost,SD and index of portfolios on the EF in the same matrix
cost_SD_idx_port_EF_temp = [cost_SD_port_EF_temp idx_port_EF_temp'];
cost_SD_idx_port_EF = -sortrows(-cost_SD_idx_port_EF_temp); 
% The above command is needed since function prtp might not return the portfolios
% on the EF in the correct order!!

% Emissions of the efficient portfolios
CO2_EF = totalCO2(cost_SD_idx_port_EF(:,3)); 
NOX_EF = totalNOX(cost_SD_idx_port_EF(:,3));
SO2_EF = totalSO2(cost_SD_idx_port_EF(:,3)); 
PM_EF = totalPM(cost_SD_idx_port_EF(:,3)); 

cost_SD_CO2_NOX_SO2_PM_idx_port_EF = [cost_SD_idx_port_EF(:,1) ...
    cost_SD_idx_port_EF(:,2) CO2_EF NOX_EF SO2_EF PM_EF ...
    cost_SD_idx_port_EF(:,3)];

% graph in per unit ($/Mwh)
cost_SD_mix_pu = [cost_mean_mix_pu SD_pu ];
[cost_SD_port_EF_pu_temp, idx_port_EF_pu_temp] = prtp(cost_SD_mix_pu);

cost_SD_idx_port_EF_pu_temp = [cost_SD_port_EF_pu_temp idx_port_EF_pu_temp'];
cost_SD_idx_port_EF_pu = -sortrows(-cost_SD_idx_port_EF_pu_temp); 

CO2_EF_pu = totalCO2(cost_SD_idx_port_EF_pu(:,3)); % CO2 emissions of the efficient portfolios
NOX_EF_pu = totalNOX(cost_SD_idx_port_EF_pu(:,3));
SO2_EF_pu = totalSO2(cost_SD_idx_port_EF_pu(:,3));
PM_EF_pu = totalPM(cost_SD_idx_port_EF_pu(:,3));

cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu = [cost_SD_idx_port_EF_pu(:,1) ...
    cost_SD_idx_port_EF_pu(:,2) CO2_EF_pu NOX_EF_pu SO2_EF_pu PM_EF_pu ...
    cost_SD_idx_port_EF_pu(:,3)];

%% average FC and VC for efficient generation portfolios
avg_FC_VC_port_EF = ...
    avg_FC_VC_port(cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu(:,7),:);
CF_mean_port_EF = ...
    CF_mean_port(cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu(:,7),:);

avg_FC_VC_percent_port = (avg_FC_VC_port./ ...
    (repmat(sum(avg_FC_VC_port,2),1,2)))*100;

%% ***********************************************************************
% Percentage share of generation technology and installed capacity
% of optimal generation portfolios on efficient portfolios
% ************************************************************************
% Share of install capacity
IC_percent_port_total_EF = ...
    IC_percent_port_total(cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu(:,7),:);
IC_percent_port_withexist_total_EF = IC_percent_port_withexist_total ...
    (cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu(:,7),:);
IC_port_total_EF_GW = ...
    Installed_cap_port_total(cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu(:,7),:); % in GW

% Share of generation (energy) 
Energy_PV_allport = Energy_actual_PV_all * ones(length(mean_energy_fossil),1);
Energy_wind_allport = Energy_actual_wind_all * ...
    ones(length(mean_energy_fossil),1);
Energy_hydro_allport = Energy_hydro_project_2030 * ...
    ones(length(mean_energy_fossil),1);

% combine mean energy of new of existing generation (TWh) by each
% technology type ('coal','CCGT','Nuke','IGCC')
Energy_mean_tech_type_fossil = ...
    [(mean_energy_fossil(:,1) + mean_energy_fossil(:,5))...
    (mean_energy_fossil(:,2) + mean_energy_fossil(:,6)) ...
    (mean_energy_fossil(:,3) + mean_energy_fossil(:,7)) ...
    mean_energy_fossil(:,4)]/10^3; % in TWh

% Mean energy for each technology for all portfolios and those on EF(MWh)
% 'coal','CCGT','Nuke','IGCC','PV','Wind','Hydro'
Energy_mean_tech_type_all = [Energy_mean_tech_type_fossil ...
    Energy_PV_allport Energy_wind_allport Energy_hydro_allport]; %(TWh)
Energy_mean_tech_type_all_EF = Energy_mean_tech_type_all ...
    (cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu(:,7),:); %(TWh)

%% ***********************************************************************
% Capacity factor of each technology (combined both new and existing)
%*************************************************************************
% CF of hydro
CF_hydro = ((Energy_hydro_project_2030 * 10^3)/(IC_hydro_2030 * 8760))*100;

CF_tech_all_port = ((Energy_mean_tech_type_all*10^3)./ ...
    (Installed_cap_port_total*8760))*100;
CF_tech_all_port(isnan(CF_tech_all_port))=0;
CF_tech_all_EF = CF_tech_all_port(cost_SD_CO2_NOX_SO2_PM_idx_port_EF_pu(:,7),:);

%% ***********************************************************************
% cost and emissions efficient frontier
% ************************************************************************
% total emissions of each generation portfolios

% Cost and CO2
cost_CO2_mix = [cost_mean_mix_billion totalCO2];
[cost_CO2_port_EF_temp, idx_port_EF_CO2_temp] = prtp(cost_CO2_mix);
cost_CO2_idx_port_EF_temp = [cost_CO2_port_EF_temp idx_port_EF_CO2_temp'];
cost_CO2_idx_port_EF = -sortrows(-cost_CO2_idx_port_EF_temp);

cost_CO2_pu_mix = [cost_mean_mix_pu totalCO2];
[cost_CO2_port_EF_pu_temp, idx_port_EF_CO2_pu_temp] = prtp(cost_CO2_pu_mix);
cost_CO2_idx_port_EF_pu_temp = [cost_CO2_port_EF_pu_temp ...
    idx_port_EF_CO2_pu_temp'];
cost_CO2_idx_port_EF_pu = -sortrows(-cost_CO2_idx_port_EF_pu_temp);

% Cost and NOX
cost_NOX_mix = [cost_mean_mix_billion totalNOX];
[cost_NOX_port_EF_temp, idx_port_EF_NOX_temp] = prtp(cost_NOX_mix);
cost_NOX_idx_port_EF_temp = [cost_NOX_port_EF_temp idx_port_EF_NOX_temp'];
cost_NOX_idx_port_EF = -sortrows(-cost_NOX_idx_port_EF_temp);

cost_NOX_pu_mix = [cost_mean_mix_pu totalNOX];
[cost_NOX_port_EF_pu_temp, idx_port_EF_NOX_pu_temp] = prtp(cost_NOX_pu_mix);
cost_NOX_idx_port_EF_pu_temp = [cost_NOX_port_EF_pu_temp ...
    idx_port_EF_NOX_pu_temp'];
cost_NOX_idx_port_EF_pu = -sortrows(-cost_NOX_idx_port_EF_pu_temp);

% Cost and SO2
cost_SO2_mix = [cost_mean_mix_billion totalSO2];
[cost_SO2_port_EF_temp, idx_port_EF_SO2_temp] = prtp(cost_SO2_mix);
cost_SO2_idx_port_EF_temp = [cost_SO2_port_EF_temp idx_port_EF_SO2_temp'];
cost_SO2_idx_port_EF = -sortrows(-cost_SO2_idx_port_EF_temp);

cost_SO2_pu_mix = [cost_mean_mix_pu totalSO2];
[cost_SO2_port_EF_pu_temp, idx_port_EF_SO2_pu_temp] = prtp(cost_SO2_pu_mix);
cost_SO2_idx_port_EF_pu_temp = [cost_SO2_port_EF_pu_temp ...
    idx_port_EF_SO2_pu_temp'];
cost_SO2_idx_port_EF_pu = -sortrows(-cost_SO2_idx_port_EF_pu_temp);

% Cost and PM
cost_PM_mix = [cost_mean_mix_billion totalPM];
[cost_PM_port_EF_temp, idx_port_EF_PM_temp] = prtp(cost_PM_mix);
cost_PM_idx_port_EF_temp = [cost_PM_port_EF_temp idx_port_EF_PM_temp'];
cost_PM_idx_port_EF = -sortrows(-cost_PM_idx_port_EF_temp);

cost_PM_pu_mix = [cost_mean_mix_pu totalPM];
[cost_PM_port_EF_pu_temp, idx_port_EF_PM_pu_temp] = prtp(cost_PM_pu_mix);
cost_PM_idx_port_EF_pu_temp = [cost_PM_port_EF_pu_temp ...
    idx_port_EF_PM_pu_temp'];
cost_PM_idx_port_EF_pu = -sortrows(-cost_PM_idx_port_EF_pu_temp);

%% ***********************************************************************
% save results to a mat file
% ************************************************************************
clear 'random_RLDC'; % to exclude this matrix from the saved file to reduce the file size
clear 'overallcost_mix_withPV_Wind_Hydro';

savename = ['D:\My Documents\Aust-China work\MATLAB\Results\scn_5\' ...
    'sum_China_',num2str(PV_wind_pen(1)),'PV_',num2str(PV_wind_pen(2)), ...
        'wind_scn5_full_cleancoal.mat'];
save(savename)
save(savename)

    end
end

