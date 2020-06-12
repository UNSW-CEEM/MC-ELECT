function random_cap_cost_all = generate_cap_cost_China(n)

% This function generate n random captical cost ($/MW) of generation
% technologies. Mean capital cost and standard deviation need to be
% specified

% Using lognormal distribution to represent uncertainties

%User input for mean and SD of carbon price

%% The mean and SD are obtained from the IEA/NEA "Projected cost 
% electricity 2010"

mean_cap_cost_coal = 510000; %overnight capital cost of coal ($/MW)
mean_cap_cost_CCGT = 391000; %overnight capital cost of CCGT ($/MW)
mean_cap_cost_nuke = 3200000; %overnight capitcatl cost of OCGT ($/MW)
mean_cap_cost_IGCC = 1206000;
mean_cap_cost_PV = 830000;
mean_cap_cost_wind = 530000;
mean_cap_cost_hydro = 982000;

SD_percent_coal = 0.2;
SD_percent_CCGT = 0.1;
SD_percent_nuke = 0.3;
SD_percent_IGCC = 0.2;
SD_percent_wind = 0.2;
SD_percent_PV = 0.2;
SD_percent_hydro = 0.1;

%% Standard deviation of capital cost of each technology
SD_cap_coal = mean_cap_cost_coal * SD_percent_coal; 
SD_cap_CCGT = mean_cap_cost_CCGT * SD_percent_CCGT; 
SD_cap_nuke = mean_cap_cost_nuke * SD_percent_nuke;
SD_cap_IGCC = mean_cap_cost_IGCC * SD_percent_IGCC;
SD_cap_PV = mean_cap_cost_PV * SD_percent_PV;
SD_cap_wind = mean_cap_cost_wind * SD_percent_wind;
SD_cap_hydro = mean_cap_cost_hydro * SD_percent_hydro;

mean_cap_cost = [mean_cap_cost_coal mean_cap_cost_CCGT ...
    mean_cap_cost_nuke mean_cap_cost_IGCC mean_cap_cost_PV ...
    mean_cap_cost_wind mean_cap_cost_hydro];
SD_cap_cost = [SD_cap_coal SD_cap_CCGT SD_cap_nuke SD_cap_IGCC ...
    SD_cap_PV SD_cap_wind SD_cap_hydro];

%% convert distribution of capital cost to  lognormal
% Coal
var_cap_coal = SD_cap_coal^2;
mu_cap_coal = log((mean_cap_cost_coal^2)/...
    sqrt(var_cap_coal + mean_cap_cost_coal^2));
sigma_cap_coal = sqrt(log(var_cap_coal/...
    (mean_cap_cost_coal^2) + 1));

% CCGT
var_cap_CCGT = SD_cap_CCGT^2;
mu_cap_CCGT = log((mean_cap_cost_CCGT^2)/sqrt(var_cap_CCGT + ...
    mean_cap_cost_CCGT^2));
sigma_cap_CCGT = sqrt(log(var_cap_CCGT/(mean_cap_cost_CCGT^2)+1));

% Nuclear
var_cap_nuke = SD_cap_nuke^2;
mu_cap_nuke = log((mean_cap_cost_nuke^2)/sqrt(var_cap_nuke + ... 
    mean_cap_cost_nuke^2));
sigma_cap_nuke = sqrt(log(var_cap_nuke/(mean_cap_cost_nuke^2)+1));

% IGCC
var_cap_IGCC = SD_cap_IGCC^2;
mu_cap_IGCC = log((mean_cap_cost_IGCC^2)/sqrt(var_cap_IGCC + ... 
    mean_cap_cost_IGCC^2));
sigma_cap_IGCC = sqrt(log(var_cap_IGCC/(mean_cap_cost_IGCC^2)+1));

% PV
var_cap_PV = SD_cap_PV^2;
mu_cap_PV = log((mean_cap_cost_PV^2)/sqrt(var_cap_PV + ...
    mean_cap_cost_PV^2));
sigma_cap_PV = sqrt(log(var_cap_PV/(mean_cap_cost_PV^2)+1));

% Wind
var_cap_wind = SD_cap_wind^2;
mu_cap_wind = log((mean_cap_cost_wind^2)/sqrt(var_cap_wind + ...
    mean_cap_cost_wind^2));
sigma_cap_wind = sqrt(log(var_cap_wind/(mean_cap_cost_wind^2)+1));

% Hydro
var_cap_hydro = SD_cap_hydro^2;
mu_cap_hydro = log((mean_cap_cost_hydro^2)/sqrt(var_cap_hydro + ...
    mean_cap_cost_hydro^2));
sigma_cap_hydro = sqrt(log(var_cap_hydro/(mean_cap_cost_hydro^2)+1));

Mu_cap = [mu_cap_coal mu_cap_CCGT mu_cap_nuke mu_cap_IGCC ...
    mu_cap_PV mu_cap_wind mu_cap_hydro];
Sigma_cap = [sigma_cap_coal sigma_cap_CCGT sigma_cap_nuke sigma_cap_IGCC...
    sigma_cap_PV sigma_cap_wind sigma_cap_hydro];

%% generate random capital cost
tic
for cc = 1:n
    sample_cap = lognrnd(Mu_cap,Sigma_cap);
    cap_cost_coal = sample_cap(:,1);
    cap_cost_CCGT = sample_cap(:,2);
    cap_cost_nuke = sample_cap(:,3);
    cap_cost_IGCC = sample_cap(:,4);
    cap_cost_PV = sample_cap(:,5);
    cap_cost_wind = sample_cap(:,6);
    cap_cost_hydro = sample_cap(:,7);
%     if cap_cost_Browncoal < 1800000
%         cap_cost_Browncoal = lognrnd(mu_cap_Browncoal,sigma_cap_Browncoal);
%     end
%     if cap_cost_Blackcoal < 1600000
%         cap_cost_Blackcoal = lognrnd(mu_cap_Blackcoal,sigma_cap_Blackcoal);
%         %cap_cost_Blackcoal = 1500000;
%     end
%     if cap_cost_CCGT < 950000
%         cap_cost_CCGT = 950000;
%     end
%     if cap_cost_OCGT < 750000
%         cap_cost_OCGT = 750000;
%     end
%     if cap_cost_PV < 2100000
%         cap_cost_PV = lognrnd(mu_cap_PV, sigma_cap_PV);
%     end
    cap_cost_coal_all(cc,:) = cap_cost_coal;    
    cap_cost_CCGT_all(cc,:) = cap_cost_CCGT;
    cap_cost_nuke_all(cc,:) = cap_cost_nuke;
    cap_cost_IGCC_all(cc,:) = cap_cost_IGCC;
    cap_cost_PV_all(cc,:) = cap_cost_PV;
    cap_cost_wind_all(cc,:) = cap_cost_wind;
    cap_cost_hydro_all(cc,:) = cap_cost_hydro;
end

random_cap_cost_all = [cap_cost_coal_all cap_cost_CCGT_all ...
    cap_cost_nuke_all cap_cost_IGCC_all cap_cost_PV_all ...
    cap_cost_wind_all cap_cost_hydro_all];

%% Save the samples to files

savename = ['D:\My Documents\Aust-China work\MATLAB\Sample_file\' ...
    'sample_cap_cost_China.mat'];
save(savename,'random_cap_cost_all')


