%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loading Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Clearing Workspace
clear

%%% Loading Hatfield Time Series
load('YaquinaTS.mat');

%%% Loading Tide Data
load('HourlyTides.mat');

%%% Loading Discharge Data
load('YaquinaRiverDischarge.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Formatting Data on Time Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finding the first and last days of the Hatfield time series
first_day = min(yaquina_all.datenum);
last_day = max(yaquina_all.datenum);

%%% Creating a daily time grid for the Hatfield data
hourlytimegrid = floor(first_day):datenum(hours(1)):ceil(last_day);
dt = 1;

%%% Creating a running mean for Hatfield salinity data
yaquina_all.datenum_rm = hourlytimegrid;
yaquina_all.datetime_rm = datetime(yaquina_all.datenum_rm, 'ConvertFrom', 'datenum');
yaquina_all.salt_rm = runmean(yaquina_all.datenum_rm, yaquina_all.datenum, yaquina_all.salt);

%%% Creating a running mean for tide data
tides.datenum_rm = hourlytimegrid;
tides.datetime_rm = datetime(tides.datenum_rm, 'ConvertFrom', 'datenum');
tides.MSL_rm = runmean(tides.datenum_rm, tides.datenum, tides.MSL);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Visualizing Time Series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [100 100 1200 800])

ax1 = subplot(311);
hold on
plot(yaquina_all.datetime, yaquina_all.salt, 'k', 'DisplayName', 'Full Resolution');
plot(yaquina_all.datetime_rm, yaquina_all.salt_rm, 'r', 'DisplayName', '1 Hour Runmean')
hold off
xlim([yaquina_all.datetime_rm(1) yaquina_all.datetime_rm(end)]);
ylabel('Salinity (psu)');
legend();

ax2 = subplot(312);
hold on
plot(tides.datetime, tides.MSL, 'k', 'DisplayName','Full Resolution');
plot(tides.datetime_rm, tides.MSL_rm, 'b', 'DisplayName', '1 Hour Runmean')
hold off
xlim([yaquina_all.datetime_rm(1) yaquina_all.datetime_rm(end)]);
ylabel('Mean Sea Level (m)');
legend();

ax3 = subplot(313);
hold on
yyaxis left
plot(yaquina_all.datetime_rm, yaquina_all.salt_rm, 'r', 'DisplayName', 'Hatfield Salinity (1 Hour Runmean)')
yyaxis right
plot(tides.datetime_rm, tides.MSL_rm, 'b', 'DisplayName', 'Mean Sea Level (1 Hour Runmean)')
hold off
xlim([yaquina_all.datetime_rm(1) yaquina_all.datetime_rm(end)]);
legend()

linkaxes([ax1 ax2 ax3], 'x');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Conducting Lagged Correlation Analysis %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,~,s2x,~,s2y,k,~] = xcovar(yaquina_all.salt_rm, tides.MSL_rm, 6);

%%% Getting Pearson Correlation Coefficient from covariance
pearson_corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);

%%% Identifying max correlation coefficient
[M,I] = max(pearson_corr_coeff);
disp('   ');
disp('Highest Correlation Coefficient: ' + string(M));

%%% Calculating the corresponding time lag
disp('Corresponding Time Lag: ' + string(k(I)*dt) + ' Hour(s)');

clear Rxy mux s2x muy s2y k Nk M I pearson_corr_coeff

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Conducting Running Lagged Correlation Analysis %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating running correlation analysis window
time_window = 6;

%%% Initializing arrays
max_corr_coeff = NaN(size(hourlytimegrid));
time_lag = NaN(size(hourlytimegrid));
fitted_time_lag = NaN(size(hourlytimegrid));

for i = 1:length(hourlytimegrid)
    
    %%% Identifying indices for running mean
    idx = hourlytimegrid > hourlytimegrid(i)-time_window/2 & hourlytimegrid < hourlytimegrid(i)+time_window/2;
    
    %%% Extracting salinity and MSL data for running mean
    salt = yaquina_all.salt_rm(idx);
    MSL = tides.MSL_rm(idx);
    
    %%% 
    if length(salt(~isnan(salt))) < 24 || length(MSL(~isnan(MSL))) < 24
        continue
    end
    
    %%% Calculating covariance
    [Rxy,~,s2x,~,s2y,k,~] = xcovar(salt, MSL, 6);
    
    %%% Getting Pearson Correlation Coefficient from covariance
    corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);
    
    %%% Identifying max correlation coefficient
    [max_corr_coeff(i),I] = max(corr_coeff);

    if isnan(max_corr_coeff(i))
        continue
    end
    
    %%% Identifying corresponding time lag
    time_lag(i) = k(I)*dt;
    
    %%% Creating variables for parabolic fitting
    idy = I-2:1:I+2;
    corr_coeff_4_fitting = corr_coeff(idy);
    time_lag_4_fitting = k(idy)*dt;
    
    %%% Fitting parabola
    p = polyfit(time_lag_4_fitting, corr_coeff_4_fitting, 2);
    
    %%% Finding max
    fitted_time_lag(i) = roots(polyder(p));
    
%     f = polyval(p, k(I)-2:0.1:k(I)+2);
%     figure()
%     plot(k, corr_coeff)
%     hold on
%     plot(k(I)-2:0.1:k(I)+2, f);
%     xline(fitted_time_lag(i))
%     hold off 
end

clear idx salt MSL Rxy mux s2x muy s2y k Nk idy corr_coeff_4_fitting time_lag_4_fitting...
    p f i I time_window

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Visualizing Correlation Analysis Results %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [100 100 1200 800])

ax1 = subplot(411);
hold on
yyaxis left
plot(yaquina_all.datetime_rm, yaquina_all.salt_rm)
ylabel('Salinity (psu)');
yyaxis right
plot(tides.datetime_rm, tides.MSL_rm)
ylabel('Mean Sea Level');
hold off

ax2 = subplot(412);
plot(datetime(hourlytimegrid, 'ConvertFrom', 'datenum'), max_corr_coeff, 'k')
ylabel('Pearson Correlation Coefficient');
title('Max Correlation Coefficient')

ax3 = subplot(413);
scatter(datetime(hourlytimegrid, 'ConvertFrom', 'datenum'), time_lag, 'k')
hold on
plot(datetime(hourlytimegrid, 'ConvertFrom', 'datenum'), fitted_time_lag, 'r');
ylabel('Time Lag (Hours)');
title('Corresponding Time Lag');

ax4 = subplot(414);
plot(riverflow.datetime, riverflow.flow);
ylabel('Yaquina River Discharge (m^3/s)');
xlim([datetime(hourlytimegrid(1), 'ConvertFrom', 'datenum') datetime(hourlytimegrid(end), 'ConvertFrom', 'datenum')]);
title('River Discharge');

linkaxes([ax1 ax2 ax3 ax4], 'x');

clear ax1 ax2 ax3 ax4
