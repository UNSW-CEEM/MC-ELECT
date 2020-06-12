%--------------------- Codes_MCS_modelling.m------------------------------
% This script file for the general portfolio modelling. This modelling
% takes into account the existing coal and gas generation capacity
% %
%------------------- created 20/04/2015 by Peerapat V.--------------------
close all; clear all; clc


% fuel prices
mean_carbon = 20;           % Expected carbon price ($/tCO2)
mean_gas = 11.65;           % gas price ($/GJ)
mean_coal = 2.1;            % coal price ($/GJ)
mean_gas_CCGT = mean_gas;   % gas price for CCGT ($/GJ)
mean_gas_OCGT = mean_gas_CCGT * 1.2; % 20% gas price uplift for OCGT
mean_fuel_Cogen = 3.77;     % fuel price for cogen plants ($/GJ)
mean_fuel_Distill = 32.31;  % fuel price for distillate ($/GJ)

% PV and wind energy penetration in percentage
PV_wind_pen = [5 10]; 

% Percentage uncertainty of carbon price, gas price and demand
Year = 2030; % Specify the planning year
SD_carbon_percent = 0.3; % Standard deviation (SD) of carbon price in % 
SD_demand_percent = 0.05; % SD of residual peak demand in %
SD_carbon = mean_carbon * SD_carbon_percent; % SD in absolute value

% Specify installed fossil fuel capacity (MW) (excldue wind, PV, Hydro)
Installed_cap =  30030;     

% Specify cost of enegy not served
cost_ENS = 12900; % unit cost of energy not served ($/MWh)

% Specify number of random sample (number of Monte Carlo runs) and number
% of hours in each interval of the LDC (RLDCs)
interval = 1; % no. of interval in each bins in the average LDC and RLDCs
LDC_sample = 1000; % no. of samples (must correspond with the sample files)

% Define minimum generation of each conventional technology (MW)
min_CCGT = 0;       % new build CCGT
min_OCGT = 0;       % new build OCGT
min_Cogen = 0;      % exiting cogeneration 
min_Coal_exist = 0; % existing coal plants
min_CCGT_exist = 0; % existing CCGT
min_OCGT_exist = 0; % existing OCGT
min_Hydro = 0;      % Hydro
min_Distill =0;     % Distillate

min_gen_with_exist = [min_Coal_exist; min_CCGT; min_OCGT; ...
    min_CCGT_exist; min_OCGT_exist; min_Cogen; min_Distill];

% Define number of generation portfolio and increment of the share of each
% technology in the portfolio. I
gen_data = [];
for a = 0:10:100
    for b = 0:10:100
        for c = 0:10:100
            total_mix = a+b+c;
            if total_mix == 100
                gen_data = [gen_data; a b c];
            end
        end
    end
end

%% ***********************************************************************
% Load data and sample files
% ************************************************************************
% Reads in the following pre-processed input data files:
% 1. Demand and RLDC for different PV and wind penetration from the MAT
%    file named "LDC_data_2030.mat"
% 2. Random fuel and carbon prices and RLDC: 
%    "sample_fuel_carbon_price_mean_91_2030.mat" and
%    "random_RLDC_5PV_10wind_2030.mat"

% Load 'sample_fuel_carbon_price_mean20_2030' data file, which is provided 
% separately. This file has to be saved into the same directory as this 
% script file
load('sample_fuel_carbon_price_mean20_2030.mat')

% Load 'LDC_data_2030.mat' data file, which is provided separately. This 
% file has to be saved into the same directory as this script file
vars = {'Demand_NEM_2030', 'Energy_hydro', 'Energy_actual_PV_diffscn', ...
    'Energy_actual_existwind_diffscn', ...
    'Energy_actual_newwind_diffscn', ...
    'IC_PV_diffpen','IC_Wind_new_diffpen',...
    'IC_Wind_existing', 'RLDC_diffscn_15SC_lesshydro'};
load(['LDC_data_2030.mat'], vars{:}) 

% Load 'random_RLDC_5PV_10wind_2030.mat', which is provided separately. This 
% file has to be saved into the same directory as this script file
vars_random_RLDC = {'random_RLDC_lesshydro_5PV_10wind'};
load(['random_RLDC_5PV_10wind_2030.mat'],vars_random_RLDC{:});
eval(['random_RLDC', '= random_RLDC_lesshydro_' num2str(PV_wind_pen(1)) 'PV_' ...
    num2str(PV_wind_pen(2)) 'wind' ';'])
clear (['random_RLDC_lesshydro_',num2str(PV_wind_pen(1)),'PV_', ...
    num2str(PV_wind_pen(2)),'wind']);

% Assigning new variables based on the level of PV and wind penetration
if PV_wind_pen == [0 0]
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(1,:);
    Energy_actual_PV = Energy_actual_PV_diffscn(1);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(1);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(1);
    IC_PV = 0;
    IC_Wind_new = 0;
elseif PV_wind_pen == [5 10]
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(2,:);
    Energy_actual_PV = Energy_actual_PV_diffscn(2);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(2);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(2);
    IC_PV = IC_PV_diffpen(1);
    IC_Wind_new = IC_Wind_new_diffpen(1);
elseif PV_wind_pen == [10 20]
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(3,:);
    Energy_actual_PV = Energy_actual_PV_diffscn(3);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(3);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(3);
    IC_PV = IC_PV_diffpen(2);
    IC_Wind_new = IC_Wind_new_diffpen(2); 
elseif PV_wind_pen == [20 30]
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(4,:);
    Energy_actual_PV = Energy_actual_PV_diffscn(4);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(4);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(4);
    IC_PV = IC_PV_diffpen(3);
    IC_Wind_new = IC_Wind_new_diffpen(3); 
elseif PV_wind_pen == [30 40]
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(5,:);
    Energy_actual_PV = Energy_actual_PV_diffscn(5);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(5);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(5);
    IC_PV = IC_PV_diffpen(4);
    IC_Wind_new = IC_Wind_new_diffpen(4); 
elseif PV_wind_pen == [40 50]
    RLDC_15SC_lesshydro = RLDC_diffscn_15SC_lesshydro(6,:);
    Energy_actual_PV = Energy_actual_PV_diffscn(6);
    Energy_actual_existwind = Energy_actual_existwind_diffscn(6);
    Energy_actual_newwind = Energy_actual_newwind_diffscn(6);
    IC_PV = IC_PV_diffpen(5);
    IC_Wind_new = IC_Wind_new_diffpen(5); 
end
        
Energy_actual_wind_all = Energy_actual_existwind + Energy_actual_newwind; %(MWh)
demand = RLDC_15SC_lesshydro;
demand_1 = demand;
hour = length(demand);

% Combined energy of PV, wind and hydro
Energy_PV_Wind_Hydro = Energy_actual_PV + Energy_actual_wind_all + ...
    Energy_hydro; % (MWh)

%% ********************************************************************** 
% Determining the installed generation capacity 
% ***********************************************************************
% Specify existing fleet in 2030 (capacity that still exists in 2030) (MW)
IC_coal_exist = 19814; 
IC_CCGT_exist = 2758; 
IC_OCGT_exist = 7415;
IC_Hydro_exist = 7654; 
IC_Cogen_exist = 586; 
IC_Distill_exist = 171;

IC_Wind_total = IC_Wind_new + IC_Wind_existing;
IC_exist = [IC_CCGT_exist IC_OCGT_exist ...
    IC_Cogen_exist IC_Distill_exist];
IC_exist_fossil = [IC_coal_exist IC_CCGT_exist IC_OCGT_exist];
IC_exist_port = repmat(IC_exist,length(gen_data),1);
IC_exist_fossil_port = repmat(IC_exist_fossil,length(gen_data),1);

% Calculate the installed capacity of conventional generation technologies
% (CCGT and OCGT) to satisfy demand x% (i.e. 0.002% of unserved energy)
Installed_cap_fossil = Installed_cap - IC_Cogen_exist - IC_Distill_exist; %

% IC with hydro, PV and Wind
Installed_cap_withhydro = Installed_cap + IC_Hydro_exist;
Installed_cap_total = Installed_cap_withhydro + IC_PV + IC_Wind_new + ...
    IC_Wind_existing;

% Capacity of each fossil-fuel technology in each portfolio
Installed_cap_fossil_port = (gen_data/100) * Installed_cap_fossil;

% real build capacity in 2030 (capacity that incurs capital cost which are
% new CCGT and OCGT)
IC_realbuild_port = Installed_cap_fossil_port - IC_exist_fossil_port;
IC_realbuild_port(IC_realbuild_port < 0)=0;

% Portfolios which exclude those with existing coal capacity is greater 
% than the actual existing
port_exclude = find(IC_realbuild_port(:,1) > 0); 
IC_realbuild_port(port_exclude,:) = [];
Installed_cap_fossil_port(port_exclude,:) = [];
IC_exist_port(port_exclude,:) = [];
IC_exist_fossil_port(port_exclude,:) = []; 

% Installed capacity of each technology categorized into new and existing
IC_port_with_exist = [Installed_cap_fossil_port(:,1) IC_realbuild_port(:,2:3) ...
    Installed_cap_fossil_port(:,2:3) - IC_realbuild_port(:,2:3) ...
    IC_exist_port(:,3:4)];
IC_port_with_exist_withPVWind_hydro = [IC_port_with_exist ...
    repmat(IC_PV,size(IC_port_with_exist,1),1) ...
    repmat(IC_Wind_total,size(IC_port_with_exist,1),1) ...
    repmat(IC_Hydro_exist,size(IC_port_with_exist,1),1)];

% Share of technology categorised into new and existing
IC_percent_port_withexist_total = ...
    (IC_port_with_exist_withPVWind_hydro/Installed_cap_total)*100;

% IC with every technology (combining existing and new)
Installed_cap_port_total = [Installed_cap_fossil_port ...
    IC_port_with_exist_withPVWind_hydro(:,6:10)];
IC_percent_port_total = ...
    (Installed_cap_port_total/Installed_cap_total)*100;

% Probability that demand is met when using the adjusted LDC 
prob_demand_met = 1-(sum(random_RLDC(:,1) > Installed_cap)/LDC_sample);
expected_peak_demand = (sum(random_RLDC(:,1)))/LDC_sample;
SD_random_peak = std(random_RLDC(:,1));
SD_random_peak_percent = SD_random_peak/expected_peak_demand; 

%% ************************************************************************
% Technical parameters of each generation technology
%*************************************************************************
% Overnight capital cost ($/MW) for new build plants
mean_cap_cost_CCGT = 1113000; 
mean_cap_cost_OCGT = 751000;
mean_cap_cost_PV = 2197000;
mean_cap_cost_wind = 1816000;

% Fixed O&M cost of each technology ($/MW/year)
FOM_coal_exist = 55651; 
FOM_CCGT = 10000; 
FOM_OCGT = 4000; 
FOM_PV = 38000; 
FOM_wind_new = 40000; 
FOM_wind_exist = 23459;
FOM_CCGT_exist = 32307; 
FOM_OCGT_exist = 17355;
FOM_Hydro = 55998; 
FOM_Cogen = 26922; 
FOM_Distill = 14000;

FOM_with_exist = [FOM_coal_exist FOM_CCGT FOM_OCGT FOM_CCGT_exist ...
    FOM_OCGT_exist FOM_Cogen FOM_Distill]; % Put FOM of all technologies into one matrix

% Plant life of each new-bhild technology (years)
plant_life_CCGT = 30; 
plant_life_OCGT = 30;
plant_life_PV = 30; 
plant_life_wind = 30;

% Variable O&M cost ($/MWh)
VOM_coal_exist = 1.34;  % existing coal plants
VOM_CCGT = 4;           % New CCGT
VOM_OCGT = 10;          % New OCGT
VOM_CCGT_exist = 1.69;  % existing CCGT
VOM_OCGT_exist = 8.76;  % existing OCGT
VOM_Hydro = 6.91;       % Hydro
VOM_Cogen = 2.07;       % Cogeneration
VOM_Distill = 10.21;    % Distillate
VOM_PV = 0;             % PV
VOM_wind_exist = 2.65;  % existing wind
VOM_wind_new = 14;      % New wind

VOM_fossil_with_exist = [VOM_coal_exist VOM_CCGT VOM_OCGT VOM_CCGT_exist...
    VOM_OCGT_exist VOM_Cogen VOM_Distill]; % Put into the same matrix

% Heat rate of each technology (GJ/MWh)
HR_coal_exist = 10.59;  % existing coal plants
HR_CCGT = 7.27;         % New CCGT
HR_OCGT = 10.29;        % New OCGT
HR_CCGT_exist = 7.78;   % existing CCGT
HR_OCGT_exist = 12.9;   % existing OCGT
HR_Cogen = 9.09;        % Cogeneration
HR_Distill = 13.22;     % Distillate

% Emission Factor of each technology (tCO2/MWh)
emission_coal_exist = 0.962;    % existing coal plants
emission_CCGT = 0.368;          % New CCGT
emission_OCGT = 0.515;          % New OCGT
emission_CCGT_exist = 0.4;      % existing CCGT
emission_OCGT_exist = 0.665;    % existing OCGT
emission_Cogen = 0.466;         % Cogeneration
emission_Distill = 0.907;       % Distillate
emission_Hydro = 0;             % Hydro

emission_with_exist = [emission_coal_exist emission_CCGT emission_OCGT ...
    emission_CCGT_exist emission_OCGT_exist emission_Cogen ...
    emission_Distill];  % Put the emissions into the same matrix

%% ************************************************************************
% Calculate annualised fixed cost of each technology (using capital
% recovery (CRF) techniques)
% The fixed cost is not part of the Monte Carlo process.
%*************************************************************************
WACC = 0.1; % Weighted Average Cost of Capital (WACC)

% Calculate the Capital Cost Recovery Factor (CRF)
CRF_CCGT = (WACC*((1+WACC)^plant_life_CCGT))/...
    (((1+WACC)^plant_life_CCGT)-1);                 
CRF_OCGT = (WACC*((1+WACC)^plant_life_OCGT))/...
    (((1+WACC)^plant_life_OCGT)-1);
CRF_PV = (WACC*((1+WACC)^plant_life_PV))/...
    (((1+WACC)^plant_life_PV)-1);
CRF_wind = (WACC*((1+WACC)^plant_life_wind))/...
    (((1+WACC)^plant_life_wind)-1);

% Calculate annualized captial cost of each tech for the reference cost
% ($/MW)
annualized_cap_cost_CCGT_ref = mean_cap_cost_CCGT * CRF_CCGT;
annualized_cap_cost_OCGT_ref = mean_cap_cost_OCGT * CRF_OCGT;
annualized_cap_cost_PV_ref = mean_cap_cost_PV * CRF_PV;
annualized_cap_cost_wind_ref = mean_cap_cost_wind * CRF_wind;

annualized_cap_cost_fossil_ref = [0 annualized_cap_cost_CCGT_ref ...
    annualized_cap_cost_OCGT_ref]; % [coal CCGT OCGT]

% Total fixed cost (FC) of PV, wind and Hydro ($)
total_FC_PV = IC_PV * (FOM_PV + annualized_cap_cost_PV_ref);
total_FC_wind_new = IC_Wind_new * (FOM_wind_new + ...
    annualized_cap_cost_wind_ref);
total_FC_wind_exist = IC_Wind_existing * FOM_wind_exist;
total_FC_wind_all = total_FC_wind_new + total_FC_wind_exist;
total_FC_hydro = IC_Hydro_exist * FOM_Hydro;
total_FC_PV_Wind_Hydro = total_FC_PV + total_FC_wind_all + total_FC_hydro;

% Capital cost and fixed O&M (FOM) for every generation portfolio ($)
total_capcost_everyport = IC_realbuild_port.* ...
    repmat(annualized_cap_cost_fossil_ref,length(IC_exist_port),1);
total_capcost_everyport = [total_capcost_everyport ...
    zeros(length(IC_exist_port),4)];
total_FOM_fossil_everyport = IC_port_with_exist.* ...
        repmat(FOM_with_exist,length(IC_exist_port),1);

% Total fixed costs (Capex + FOM) of every generation portfolio considered
% ($)
total_FC_fossil_everyport = total_capcost_everyport + ...
    total_FOM_fossil_everyport;

%% *********************************************************************** 
% Calculate short run marginal cost (SRMC) of each generation technology.
% SRMC consists of fuel costs, variable O&M anc carbon costs. Note that
% the SRMC is determined for every sample of fuel and carbon prices (10,000
% simulations)
%*************************************************************************
% Extract 10,000 simulated fuel and carbon price from 
% "sample_fuel_carbon_price" matrix ($/GJ)
coal_price = sample_fuel_carbon_price(:,1); % ($/GJ)
gas_CCGT_price = sample_fuel_carbon_price(:,2);
gas_OCGT_price = sample_fuel_carbon_price(:,3);
carbon_price = sample_fuel_carbon_price(:,4);

% Calculate fuel cost of each technolgy for every sample ($/MWh)
fuel_cost_coal_exist = coal_price * HR_coal_exist;
fuel_cost_CCGT = gas_CCGT_price * HR_CCGT;
fuel_cost_OCGT = gas_OCGT_price * HR_OCGT;
fuel_cost_CCGT_exist = gas_CCGT_price * HR_CCGT_exist;
fuel_cost_OCGT_exist = gas_OCGT_price * HR_OCGT_exist;
fuel_cost_Cogen = ones(LDC_sample,1) * (mean_fuel_Cogen * HR_Cogen);
fuel_cost_Distill = ones(LDC_sample,1) * (mean_fuel_Distill * HR_Distill);
fuel_cost_Hydro = zeros(LDC_sample,1);

fuel_cost_with_exist = [fuel_cost_coal_exist fuel_cost_CCGT ...
    fuel_cost_OCGT fuel_cost_CCGT_exist fuel_cost_OCGT_exist ...
    fuel_cost_Cogen fuel_cost_Distill]; % Put fuel costs in one matrix

% Calculate carbon cost of each technology for every sample ($/MWh)
carbon_cost_coal_exist = carbon_price * emission_coal_exist;
carbon_cost_CCGT = carbon_price * emission_CCGT;
carbon_cost_OCGT = carbon_price * emission_OCGT;
carbon_cost_CCGT_exist = carbon_price * emission_CCGT_exist;
carbon_cost_OCGT_exist = carbon_price * emission_OCGT_exist;
carbon_cost_Cogen = carbon_price * emission_Cogen;
carbon_cost_Distill = carbon_price * emission_Distill;
carbon_cost_Hydro = zeros(LDC_sample,1);

carbon_cost_with_exist = [carbon_cost_coal_exist carbon_cost_CCGT ...
    carbon_cost_OCGT carbon_cost_CCGT_exist carbon_cost_OCGT_exist ...
    carbon_cost_Cogen carbon_cost_Distill]; % Put carbon costs in one matrix

% calculate variable operating cost per unit ($/MWh) of each technology 
VOM_all_fossil_with_exist = ones(LDC_sample,1) * VOM_fossil_with_exist;

% SRMC of the fossil fuel technologies for every set of random fuel and
% carbon price
VC_pu_with_exist = VOM_all_fossil_with_exist + fuel_cost_with_exist + ...
    carbon_cost_with_exist;

% Total short run variable costs for PV, wind and Hydro ($)
total_VC_PV = (Energy_actual_PV * VOM_PV); 
total_VC_wind_new = Energy_actual_newwind * VOM_wind_new;
total_VC_wind_exist = Energy_actual_existwind * VOM_wind_exist;
total_VC_wind_all = total_VC_wind_exist + total_VC_wind_new;
total_VC_hydro = Energy_hydro * VOM_Hydro;

total_VC_PV_Wind_Hydro = total_VC_PV + total_VC_wind_all + total_VC_hydro;
total_VC_PV_Wind_Hydro_everyport = repmat(total_VC_PV_Wind_Hydro, ...
    length(IC_port_with_exist),1);
total_VC_PV_Wind_Hydro_everysample = ...
    repmat(total_VC_PV_Wind_Hydro,LDC_sample,1);

%%************************************************************************* 
% Determine merit order dispatch of each technology based on their variable
% costs. 
% *************************************************************************
%       Note that PV, wind and hydro are considered exogeneous and 
%       therefore not included in the dispatch (assume priority dispatch 
%       of those technologies based on historical generation data)

% Technologies that are subject to dispatch
tech = {'existing coal','CCGT','OCGT','existing CCGT','existing OCGT', ...
    'Cogen', 'Distillate'};

% Loop through to determine dispatch for every generation portfolio
% considered, and for every set of random fuel and carbon prices. 
for r = 1:length(IC_exist_port)                     % For every generation portfolio 
    IC_with_exist = IC_port_with_exist(r,:);
    IC_realbuild = IC_realbuild_port(r,:);
    cap_tech = IC_with_exist - min_gen_with_exist'; % available capacity of each technology for the dispatch
    gen_output = zeros(length(cap_tech),length(demand)); % initialize matrix
    total_FC_fossil = repmat(total_FC_fossil_everyport(r,:),LDC_sample,1);
    total_FC = [total_FC_fossil ...                 % total fixed cost of all the technologies
        repmat(total_FC_PV_Wind_Hydro,LDC_sample,1)]; 
    for k = 1:size(random_RLDC,1)                   % For each set of random parameters (fuel, carbon prices and LDC)
        LDC = random_RLDC(k,:);
        LDC_1 = LDC;
        [VC_sort,idx] = sort(VC_pu_with_exist(k,:)); % Rank each technology based on their SRMCs
        cap_sort = cap_tech(idx);                   % capacity of each tech after sorted based on VC
        cap(k,:) = cap_sort;
        for i=1:length(cap_tech)                    % Determine output(MW) of each technology in each segment of the LDC
            gen_output(i,:) = min(LDC,cap_sort(i).*ones(1,hour));
            LDC = LDC - gen_output(i,:);
            [sortidx, ind] = sort(idx);             % sort the index to arrange back to the original order
            gen_tech = gen_output(ind,:);           % show the generation output in the order coal->CCGT->OCGT
        end
        total_gen = sum(gen_tech,1);                % sum of generation output from every technology
        ENS_period = sum((total_gen*1.000000000009) < LDC_1);       % Determine the number of interval that load is not met
        ENS_period_all(k,:) = ENS_period;
        ENS = (LDC_1(1:ENS_period) - Installed_cap)*interval;       % amount of energy not served in a year (MWh) which is the same for every portfolio
        ENS_all(k,:) = sum(ENS);
        ENS_cost = sum(ENS * cost_ENS);                             % total cost of unserved energy ($)
        ENS_cost_all(k,:) = (ENS_cost/1000000);                     % ($Million)
        energy_tech_fossil = gen_tech * interval;                   % energy of each tech in each segment of each  LDC
        energy_tech_year_fossil = sum(energy_tech_fossil,2);        % energy of each technology in a year for a set of random parameters
        energy_fossil(k,:) = energy_tech_year_fossil;               % annual energy of each tech. for every sample of fuel and carbon price
        VC_fossil = energy_tech_year_fossil'.*VC_pu_with_exist(k,:);% VC of each tech for each carbon price
        VC_tech_fossil(k,:) = VC_fossil;
        CF = (energy_tech_year_fossil'./(IC_with_exist * length(LDC_1) ...
            * interval)) * 100;
        CF_all(k,:) = CF;                                           % Capacity factor (CF) of each tech for every set of random parameters
        CF_all(isnan(CF_all))=0;
        CF_mean = mean(CF_all,1);                                   % Mean CF for a generation portfolio
        LDC = LDC_1;
    end
    CF_mean_port(r,:) = CF_mean;
    VC_tech_withPV_Wind_Hydro = [VC_tech_fossil ...                 % VC for every set of random parameters
          total_VC_PV_Wind_Hydro_everysample];                         
    VC_mean_port_fossil(r,:) = sum(VC_tech_fossil,1)/LDC_sample;    % mean VC by each technology for every generation portfolio
    FC_mean_port_fossil(r,:) = sum(total_FC_fossil,1)/LDC_sample;
    totalcost_fossil = VC_tech_fossil + total_FC_fossil;
    totalcost_withPV_Wind_Hydro = VC_tech_withPV_Wind_Hydro + total_FC; 
    totalcost_tech_fossil = totalcost_fossil/1000000;               % total cost of each technology 
    totalcost_tech_withPV_Wind_Hydro = totalcost_withPV_Wind_Hydro/1000000; % $ Million
    mean_energy_fossil(r,:) = sum(energy_fossil,1)/LDC_sample;              % Mean energy for every generation portfolio
    CO2_tech = (repmat(emission_with_exist,[LDC_sample,1])).* energy_fossil;
    total_CO2 = sum(CO2_tech,2);                                        % total CO2 emission of generation portfolio in a year for every set of random fuel and carbon prices
    overallcost_mix_fossil(r,:) = sum(totalcost_tech_fossil,2) + ...    % overall annual generation cost of fossil fuel technologies ($Million) for every portfolio ...
        ENS_cost_all;                                                   % and every set of random fuel and carbon prices (including ENS cost)
    overallcost_mix_withPV_Wind_Hydro(r,:) = ...                         
        sum(totalcost_tech_withPV_Wind_Hydro,2) + ENS_cost_all;         % overall annual generation cost of all of the technologies ($Million) for every portfolio
end

overallcost_mix_billion = overallcost_mix_withPV_Wind_Hydro/1000;       % overall annual generation cost of all of the technologies ($Billion) for every portfolio

% Calculate CO2 emission of each generation portfolio
for u=1:length(mean_energy_fossil)
    CO2(u,:) = emission_with_exist.*mean_energy_fossil(u,:);
end
totalCO2 = sum(CO2,2)/1000000; % total CO2 emissions (M.ton) for each generation mix

% Average total variable cost and fixed cost of each technology over 10,000 
% simulationsfor every generation portfolio ($ Million)
VC_mean_port = [VC_mean_port_fossil total_VC_PV_Wind_Hydro_everyport(1:r)];
avg_VC_port = sum(VC_mean_port,2)/1000000;      % Average variable cost ($ Million)
avg_FC_PV_Wind_Hydro = total_FC_PV_Wind_Hydro;  % Average fixed cost is the same for every generation portfolio
FC_port = [FC_mean_port_fossil (avg_FC_PV_Wind_Hydro.*ones(r,1))]; 
avg_FC_port = sum(FC_port,2)/1000000;            % Average fixed costs($ Million)

avg_FC_VC_port = [avg_FC_port avg_VC_port]; % put FC and VC into one matrix
avg_totalcost = sum(avg_FC_VC_port,2); 

% Mean energy of each technology for each generation portfolio 
mean_energy_withPV_Wind_Hydro = [mean_energy_fossil ...
    (Energy_PV_Wind_Hydro * ones(size(mean_energy_fossil,1),1))];
mean_energy_total = sum(mean_energy_withPV_Wind_Hydro,2);

% Total energy demand (MWh) to be met by fossil-fuel gen and PV in a year
total_energy_fossil = sum(demand) * interval;
total_energy = total_energy_fossil + Energy_PV_Wind_Hydro;

% Expected cost and SD of every generation portfolio 
cost_mean_mix = sum(overallcost_mix_withPV_Wind_Hydro,2)/LDC_sample; % ($Million)
SD_mix = std(overallcost_mix_withPV_Wind_Hydro,0,2); % ($Million)
cost_mean_mix_billion = cost_mean_mix/1000; % ($Billion)
SD_mix_billion = SD_mix/1000;               % ($Billion)

% proportion of the cost of unserved energy to total cost
mean_ENS_cost = sum(ENS_cost_all)/length(random_RLDC); % Mean cost of energy not served which is the same for every generation portfolio
mean_ENS_cost_all = mean_ENS_cost*ones(r,1);
proportion_ENS_cost = ((mean_ENS_cost_all/1000000)./cost_mean_mix)*100; %Percentage of ENS cost to the overall costs

% Energy by technology type
Energy_PV_allport = Energy_actual_PV * ones(length(mean_energy_fossil),1);
Energy_wind_allport = Energy_actual_wind_all * ...
    ones(length(mean_energy_fossil),1);
Energy_hydro_allport = Energy_hydro * ones(length(mean_energy_fossil),1);

% combine mean energy of new of existing generation (MWh)
% 'coal','CCGT','OCGT','Cogen','Distillate'
Energy_mean_tech_type_fossil = [mean_energy_fossil(:,1) ...
    (mean_energy_fossil(:,2) + mean_energy_fossil(:,4)) ...
    (mean_energy_fossil(:,3) + mean_energy_fossil(:,5)) ...
    mean_energy_fossil(:,6) mean_energy_fossil(:,7)];

% Mean energy for each technology for all portfolios 
% 'coal','CCGT','OCGT','Cogen','Distillate','PV','Wind','Hydro'
Energy_mean_tech_type_all = [Energy_mean_tech_type_fossil ...
    Energy_PV_allport Energy_wind_allport Energy_hydro_allport]/10^6; %(TWh)

% Mean capacity factor of each technology
CF_tech_all_port = ((Energy_mean_tech_type_all*10^6)./ ...
    (Installed_cap_port_total*8760))*100;
CF_tech_all_port(isnan(CF_tech_all_port))=0;

%% ***********************************************************************
% Expected (mean) annual generation cost in per unit ($/MWh) and the 
% standard deviation (SD) of each generation portfolio
%*************************************************************************
energy_random_RLDC = sum(random_RLDC,2) * interval;
energy_all_random_RLDC = energy_random_RLDC + Energy_actual_PV + ...
    Energy_actual_wind_all + Energy_hydro;
overallcost_mix_pu = (overallcost_mix_withPV_Wind_Hydro*10^6)./ ...
    repmat(energy_all_random_RLDC,1,length(IC_exist_port))';
cost_mean_mix_pu = sum(overallcost_mix_pu,2)/LDC_sample; % ($/MWh)
SD_pu = std(overallcost_mix_pu,0,2);

%% *********************************************************************** 
% Identify generation portfolios on the cost-risk efficient frontier (EF)
% by calling the function "prtp"  to identify portfolios along the efficient
% Function "prtp" is provided as a separate file and must be saved in the
% same directory as this scrpit file
% ************************************************************************
cost_SD_mix_pu = [cost_mean_mix_pu SD_pu ];
[cost_SD_port_EF_pu_temp, idx_port_EF_pu_temp] = prtp(cost_SD_mix_pu);
cost_SD_idx_port_EF_pu_temp = [cost_SD_port_EF_pu_temp idx_port_EF_pu_temp'];
cost_SD_idx_port_EF_pu = -sortrows(-cost_SD_idx_port_EF_pu_temp);
CO2_EF_pu = totalCO2(cost_SD_idx_port_EF_pu(:,3));             % CO2 emissions of the efficient portfolios

cost_SD_CO2_idx_port_EF_pu = [cost_SD_idx_port_EF_pu(:,1) ...
    cost_SD_idx_port_EF_pu(:,2) CO2_EF_pu cost_SD_idx_port_EF_pu(:,3)];  % put the cost, SD, CO2 and portfolio index into the same matrix.

% average FC and VC of generation portfolios on the EF
avg_FC_VC_port_EF = avg_FC_VC_port(cost_SD_CO2_idx_port_EF_pu(:,4),:);
CF_mean_port_EF = CF_mean_port(cost_SD_CO2_idx_port_EF_pu(:,4),:);
avg_FC_VC_percent_port = (avg_FC_VC_port./ ...
    (repmat(sum(avg_FC_VC_port,2),1,2)))*100;

% Percentage share of generation technology and installed capacity
% of optimal generation portfolios on efficient portfolios
IC_percent_port_total_EF = ...
    IC_percent_port_total(cost_SD_CO2_idx_port_EF_pu(:,4),:);
IC_percent_port_withexist_total_EF = ...
    IC_percent_port_withexist_total(cost_SD_CO2_idx_port_EF_pu(:,4),:);
IC_port_total_EF_GW = ...
    Installed_cap_port_total(cost_SD_CO2_idx_port_EF_pu(:,4),:)/1000; % in GW

% Mean annual energy of each technology for portfolios on the Efficient 
% Frontier 
Energy_mean_tech_type_all_EF = ...
    Energy_mean_tech_type_all(cost_SD_CO2_idx_port_EF_pu(:,4),:); %(TWh)

% Mean capacity factor for all technology (combined both new and existing
% for the same technology)
CF_tech_all_EF = CF_tech_all_port(cost_SD_CO2_idx_port_EF_pu(:,4),:);

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
port_cost_SD_CO2_idx = [cost_mean_mix_billion(port_cost_index) ...
    SD_mix_billion(port_cost_index) totalCO2(port_cost_index) ...
    port_cost_index];
port_cost_SD_CO2_pu_idx = [cost_mean_mix_pu(port_cost_pu_index) ...
    SD_pu(port_cost_pu_index) totalCO2(port_cost_pu_index) ...
    port_cost_pu_index];

% matrix containing cost, SD, CO2 of the portfolio
cost_SD_CO2_idx_port = [cost_mean_mix_billion SD_mix_billion ...
    totalCO2];
cost_SD_CO2_idx_pu_port = [cost_mean_mix_pu SD_pu totalCO2];
  
%% ************************************************************************ 
% save results to a mat file
%*************************************************************************

% Save the following matrix to produce output results
% 'PV_wind_pen' - PV and wind energy penetration in (%)
% 'Demand_NEM_2030' - hourly demand in 2030 (MW)
% 'RLDC_15SC_lesshydro' - Residual load duration cuve after accounted for
%                         15% minimum synchronous constraint and hydro
% 'overallcost_mix_pu' - annual generation cost of every portfolio for each
%                        of the 10,000 simulations
% 'cost_SD_CO2_idx_pu_port' - A result matrix containing expected costs, 
%                             standard deviation (SD) of costs, 
%                             mean CO2 emissions of every generation
%                             portfolio considered in the simulation
% 'IC_percent_port_withexist_total' - percentage of capacity share of each 
%                                     technology in each of the portfolios 
% 'IC_port_with_exist_withPVWind_hydro' - installed capacity of each
%                                         technology in each of the 
%                                         portfolios considered
% 'cost_SD_CO2_idx_port_EF_pu' - A result matrix containing expected costs, 
%                                standard deviation (SD) of costs, 
%                                mean CO2 emissions of generation
%                                portfolio on the Efficient Frontier
% 'IC_percent_port_total_EF' - percentage of capacity share of each 
%                              technology in the portfolios on the Frontier
% 'IC_percent_port_withexist_total_EF' - percentage of capacity share of each 
%                                        technology (separated into new 
%                                        and existing) in the portfolios on
%                                        the Frontier
% 'IC_port_total_EF_GW' - installed capacity of eacn technology in the 
%                         portfolios on the Efficient Frontier
% 'Energy_mean_tech_type_all_EF' - mean energy of each technology in the
%                                  portfolios on the Efficient Frontier
% 'CF_tech_all_EF' - mean capacity factor of each technology in the
%                    portfolio on the Efficient Frontier

result_vars = {'PV_wind_pen', 'Demand_NEM_2030', 'RLDC_15SC_lesshydro', ...
    'overallcost_mix_pu','cost_SD_CO2_idx_pu_port', ...
    'IC_percent_port_withexist_total', ...
    'IC_port_with_exist_withPVWind_hydro', ...
    'cost_SD_CO2_idx_port_EF_pu',  ...
    'IC_percent_port_total_EF', 'IC_percent_port_withexist_total_EF', ...
    'IC_port_total_EF_GW', 'Energy_mean_tech_type_all_EF', ...
    'CF_tech_all_EF'};
 
savename = ['D:\test.mat']; % Users can change the directory and file name 
save(savename,result_vars{:})




