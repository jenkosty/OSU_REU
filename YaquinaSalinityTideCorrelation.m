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
yaquina_all.datenum_rm1hr = hourlytimegrid;
yaquina_all.datetime_rm1hr = datetime(yaquina_all.datenum_rm1hr, 'ConvertFrom', 'datenum');
yaquina_all.salt_rm1hr = runmean(yaquina_all.datenum_rm1hr, yaquina_all.datenum, yaquina_all.salt);

%%% Creating a running mean for tide data
tides.datenum_rm1hr = hourlytimegrid;
tides.datetime_rm1hr = datetime(tides.datenum_rm1hr, 'ConvertFrom', 'datenum');
tides.MSL_rm1hr = runmean(tides.datenum_rm1hr, tides.datenum, tides.MSL);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Visualizing Time Series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [100 100 1200 800])

ax1 = subplot(311);
hold on
plot(yaquina_all.datetime, yaquina_all.salt, 'k', 'DisplayName', 'Full Resolution');
plot(yaquina_all.datetime_rm1hr, yaquina_all.salt_rm1hr, 'r', 'DisplayName', '1 Hour Runmean')
hold off
xlim([yaquina_all.datetime_rm1hr(1) yaquina_all.datetime_rm1hr(end)]);
ylabel('Salinity (psu)');
legend();

ax2 = subplot(312);
hold on
plot(tides.datetime, tides.MSL, 'k', 'DisplayName','Full Resolution');
plot(tides.datetime_rm1hr, tides.MSL_rm1hr, 'b', 'DisplayName', '1 Hour Runmean')
hold off
xlim([yaquina_all.datetime_rm1hr(1) yaquina_all.datetime_rm1hr(end)]);
ylabel('Mean Sea Level (m)');
legend();

ax3 = subplot(313);
hold on
yyaxis left
plot(yaquina_all.datetime_rm1hr, yaquina_all.salt_rm1hr, 'r', 'DisplayName', 'Hatfield Salinity (1 Hour Runmean)')
yyaxis right
plot(tides.datetime_rm1hr, tides.MSL_rm1hr, 'b', 'DisplayName', 'Mean Sea Level (1 Hour Runmean)')
hold off
xlim([yaquina_all.datetime_rm1hr(1) yaquina_all.datetime_rm1hr(end)]);
legend()

linkaxes([ax1 ax2 ax3], 'x');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Conducting Lagged Correlation Analysis %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,mux,s2x,muy,s2y,k,Nk] = xcovar(yaquina_all.salt_rm1hr, tides.MSL_rm1hr, 6);

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

for i = 1:length(hourlytimegrid)
    idx = hourlytimegrid > hourlytimegrid(i)-time_window/2 & hourlytimegrid < hourlytimegrid(i)+time_window/2;
    time = hourlytimegrid(idx);
    salt = yaquina_all.salt_rm1hr(idx);
    MSL = tides.MSL_rm1hr(idx);
    
    %%% Calculating covariance
    [Rxy,mux,s2x,muy,s2y,k,Nk] = xcovar(salt, MSL, 6);
    
    %%% Getting Pearson Correlation Coefficient from covariance
    pearson_corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);
    
    %%% Identifying max correlation coefficient
    [max_corr_coeff(i),I] = max(pearson_corr_coeff);
    
    %%% Identifying corresponding time lag
    time_lag(i) = k(I)*dt;
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Renderer', 'painters', 'Position', [100 100 1200 800])

ax1 = subplot(411);
hold on
yyaxis left
plot(yaquina_all.datetime_rm1hr, yaquina_all.salt_rm1hr)
ylabel('Salinity (psu)');
yyaxis right
plot(tides.datetime_rm1hr, tides.MSL_rm1hr)
ylabel('Mean Sea Level');
hold off

ax2 = subplot(412);
plot(datetime(hourlytimegrid, 'ConvertFrom', 'datenum'), max_corr_coeff, 'k')
ylabel('Pearson Correlation Coefficient');
title('Max Correlation Coefficient')

ax3 = subplot(413);
scatter(datetime(hourlytimegrid, 'ConvertFrom', 'datenum'), time_lag, 'k')
ylabel('Time Lag (Hours)');
title('Corresponding Time Lag');

ax4 = subplot(414);
plot(riverflow.datetime, riverflow.flow);
ylabel('Yaquina River Discharge (m^3/s)');
xlim([datetime(hourlytimegrid(1), 'ConvertFrom', 'datenum') datetime(hourlytimegrid(end), 'ConvertFrom', 'datenum')]);
title('River Discharge');

linkaxes([ax1 ax2 ax3 ax4], 'x');
