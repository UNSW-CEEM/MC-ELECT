function sample_fuel_carbon_price = generate_coal_gas_nuke_carbon_China_1(LDC_sample)

% This function generate n samples of correlated coal, gas, nuclear fuel 
% and carbon prices. Mean carbon price and standard deviation need to be 
% specified

% User input mean carbon price
mean_carbon = input('Input mean carbon price ($/ton of CO2) = '); 
mean_gas = input('Input mean gas price ($/GJ) = ');
mean_coal = input('Input mean coal price ($/GJ) = ');

SD_carbon_percent = 0.5; % Percentage uncertainty of carbon price
mean_carbon(mean_carbon == 0) = 0.0000000001; % to prevent the case of zero carbon price
SD_carbon = mean_carbon * SD_carbon_percent; 

%mean_gas = input('Input mean gas price ($/GJ) =');
n = LDC_sample;

%% Fuel prices ($/GJ) and their SD as a percentage of the mean
% Mean fuel costs
% mean_coal = 4.52;
% mean_gas = 8.938;
mean_nuke = 0.7;

% Standard deviation as % of mean costs
SD_coal_percent = 0.1;
SD_gas_percent = 0.3;
SD_nuke_percent = 0.05;

%% Define no.of samples, mean, SD of coal, gas, carbon prices and 
% correlations

% Mean and SD of fuel prices ($/GJ)
SD_coal = SD_coal_percent * mean_coal; 
SD_gas = SD_gas_percent * mean_gas; 
SD_nuke = SD_nuke_percent * mean_nuke;

SD = [SD_coal; SD_gas; SD_nuke; SD_carbon];
mean = [mean_coal; mean_gas; mean_nuke; mean_carbon];

% Define correlation factor between gas, coal, carbon prices
corr_coal_gas = 0.6; % correlation between coal and gas price
corr_coal_nuke = 0; % correlation between coal and nuclear fuel price
corr_gas_nuke = 0; % correlation between gas and nuclear fuel price
corr_gas_carbon = 0.45; % correlation between gas and carbon price
corr_coal_carbon = -0.35; % correlation between coal and carbon price
corr_nuke_carbon = 0;


% Establish Correlation matrix
corr = [1 corr_coal_gas corr_coal_nuke corr_coal_carbon; ...
    corr_coal_gas 1 corr_gas_nuke corr_gas_carbon;...
    corr_coal_nuke corr_gas_nuke 1 corr_nuke_carbon; ...
    corr_coal_carbon corr_gas_carbon corr_nuke_carbon 1];

% Determine covariance between uncertain variables
cov = corr2cov(SD,corr); 

%% convert to distribution to lognormal
var_carbon = SD_carbon^2;
mu_carbon = log((mean_carbon^2)/sqrt(var_carbon + mean_carbon^2));
sigma_carbon = sqrt(log(var_carbon/(mean_carbon^2) + 1));

%convert gas price distribution to lognormal
var_gas = SD_gas^2;
mu_gas = log((mean_gas^2)/sqrt(var_gas + mean_gas^2));
sigma_gas = sqrt(log(var_gas/(mean_gas^2) + 1));

%convert coal price distribution to lognormal
var_coal = SD_coal^2;
mu_coal = log((mean_coal^2)/sqrt(var_coal +...
    mean_coal^2));
sigma_coal = sqrt(log(var_coal/(mean_coal^2) + 1));

%convert nuclear fuel price distribution to lognormal
var_nuke = SD_nuke^2;
mu_nuke = log((mean_nuke^2)/sqrt(var_nuke +...
    mean_nuke^2));
sigma_nuke = sqrt(log(var_nuke/(mean_nuke^2) + 1));

Mu = [mu_coal mu_gas mu_nuke mu_carbon];
Sigma = [sigma_coal sigma_gas sigma_nuke sigma_carbon];


%% Call function MVLogNRand to generate correlated random gas, coal and 
%carbon prices from theri multivariate lognormal distribution

sample_fuel_carbon_price = MvLogNRand(Mu,Sigma,n,corr); 

%y = sample_fuel_carbon_price;

savename = ['D:\My Documents\Aust-China work\MATLAB\Sample_file\' ...
    'sample_fuel_carbon_price_mean_',num2str(mean_carbon),'_gasprice_', ...
    num2str(mean_gas),'_coalprice_',num2str(mean_coal),'.mat'];
save(savename,'sample_fuel_carbon_price')

