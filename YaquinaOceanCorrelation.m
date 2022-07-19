%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loading Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Clearing Workspace
clear

%%% Loading Hatfield Time Series
load('YaquinaTS.mat');

%%% Loading OISM Time Series
load('OOI_CE_OISM.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Finding Max Daily Salinity in Hatfield %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finding the first and last days of the Hatfield time series
first_day = min(yaquina_all.datenum);
last_day = max(yaquina_all.datenum);

%%% Creating a daily time grid for the Hatfield data
dailytimegrid = floor(first_day):1:ceil(last_day);

clear first_day last_day

%%% Finding the max daily salinity
for i = length(dailytimegrid)-1:-1:1
    
    %%% Daily indices
    idx = yaquina_all.datenum > dailytimegrid(i) & yaquina_all.datenum < dailytimegrid(i+1);
    
    %%% Finding max daily salinity
    [hatfield_HT.salt(i), idy] = max(yaquina_all.salt(idx));
    
    %%% Finding associated temperature
    temp_idx = yaquina_all.temp(idx);
    hatfield_HT.temp(i) = temp_idx(idy);
    
    %%% Finding associated time
    times_idx = yaquina_all.datenum(idx);
    hatfield_HT.datenum(i) = times_idx(idy);
    hatfield_HT.datetime(i) = datetime(hatfield_HT.datenum(i), 'ConvertFrom', 'datenum');
end

clear dailydatetimes idx idy i temp_idx times_idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Interpolating Hatfield High Tide Data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noontimegrid = dailytimegrid+0.5;

hatfield_HT.datenum_interp = noontimegrid;
hatfield_HT.datetime_interp = datetime(hatfield_HT.datenum_interp, 'ConvertFrom', 'datenum');

%%% Interpolating salinity data
hatfield_HT.salt_interp = interp1(hatfield_HT.datenum, hatfield_HT.salt, hatfield_HT.datenum_interp); 

%%% Interpolating temperature data
hatfield_HT.temp_interp = interp1(hatfield_HT.datenum, hatfield_HT.temp, hatfield_HT.datenum_interp); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Calculating a Running Mean of OISM Data on Time Grid %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oism.datenum_1dayrunmean = noontimegrid;
oism.datetime_1dayrunmean = datetime(noontimegrid, 'ConvertFrom', 'datenum');

%%% 1m runmean for salinity and temperature
oism.salt1m_1dayrunmean = runmean(noontimegrid, datenum(oism.salt1m.datetime), oism.salt1m.salt);
oism.temp1m_1dayrunmean = runmean(noontimegrid, datenum(oism.temp1m.datetime), oism.temp1m.temp);

%%% 7m runmean for salinity and temperature
oism.salt7m_1dayrunmean = runmean(noontimegrid, datenum(oism.salt7m.datetime), oism.salt7m.salt);
oism.temp7m_1dayrunmean = runmean(noontimegrid, datenum(oism.temp7m.datetime), oism.temp7m.temp);

%%% Averaging 1m and 7m runmeans for salinity and temperature
oism.avgsalt_1dayrunmean = mean([oism.salt1m_1dayrunmean; oism.salt7m_1dayrunmean]);
oism.avgtemp_1dayrunmean = mean([oism.temp1m_1dayrunmean; oism.temp7m_1dayrunmean]);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Visualizing Time Series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Hatfield Figure

figure('Renderer', 'painters', 'Position', [100 100 1200 800])
sgtitle('Hatfield Summary Figure');

plot(yaquina_all.datetime, yaquina_all.salt, 'Color', [.7 .7 .7], 'DisplayName', 'Hatfield Full Resolution Salinity');
hold on
plot(vertcat(hatfield_HT.datetime), vertcat(hatfield_HT.salt), 'k', 'DisplayName', 'Max Daily Salinity');
ylabel('Salinity (psu)');
legend();

%%% OISM Figure

figure('Renderer', 'painters', 'Position', [100 100 1200 800])
sgtitle('OISM Summary Figure');

ax1 = subplot(311);
plot(oism.salt1m.datetime, oism.salt1m.salt, 'Color', [.7 .7 .7], 'DisplayName','OISM Full Resolution (1m)')
hold on
plot(oism.datetime_1dayrunmean, oism.salt1m_1dayrunmean, 'r', 'DisplayName', 'OISM 1 Day Runmean (1m)');
xlim([oism.datetime_1dayrunmean(1) oism.datetime_1dayrunmean(end)])
ylabel('Salinity (psu)');
legend();

ax2 = subplot(312);
plot(oism.salt7m.datetime, oism.salt7m.salt, 'Color', [.7 .7 .7], 'DisplayName','OISM Full Resolution (7m)')
hold on
plot(oism.datetime_1dayrunmean, oism.salt7m_1dayrunmean, 'b', 'DisplayName', 'OISM 1 Day Runmean (7m)');
xlim([oism.datetime_1dayrunmean(1) oism.datetime_1dayrunmean(end)])
ylabel('Salinity (psu)');
legend();

ax3 = subplot(313);
plot(oism.datetime_1dayrunmean, oism.salt1m_1dayrunmean, 'r', 'DisplayName', 'OISM 1 Day Runmean (1m)');
hold on
plot(oism.datetime_1dayrunmean, oism.salt7m_1dayrunmean, 'b', 'DisplayName', 'OISM 1 Day Runmean (7m)');
plot(oism.datetime_1dayrunmean, oism.avgsalt_1dayrunmean, 'k', 'DisplayName', 'OISM 1 Day Runmean (Average)');
xlim([oism.datetime_1dayrunmean(1) oism.datetime_1dayrunmean(end)])
ylabel('Salinity (psu)');
legend()

linkaxes([ax1 ax2 ax3], 'x');

%%% OISM / Hatfield Figure

figure('Renderer', 'painters', 'Position', [100 100 1200 800])
sgtitle('Hatfield and OISM Correlation Summary Figure')

ax1 = subplot(211);
plot(hatfield_HT.datetime_interp, hatfield_HT.salt_interp, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Hatfield High Tide')
hold on
plot(oism.datetime_1dayrunmean,oism.avgsalt_1dayrunmean, 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'OISM 1m/7m Average')
ylabel('Salinity (psu)');
legend()

ax2 = subplot(212);
plot(hatfield_HT.datetime_interp, hatfield_HT.temp_interp, 'Color', [0 0.4470 0.7410], 'DisplayName', 'Hatfield High Tide')
hold on
plot(oism.datetime_1dayrunmean,oism.avgtemp_1dayrunmean, 'Color', [0.9290 0.6940 0.1250], 'DisplayName', 'OISM 1m/7m Average')
ylabel('Temperature (degC)');
legend()

linkaxes([ax1 ax2], 'x')

clear ax1 ax2 ax3
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Conducting Lagged Correlation Analysis %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dry_start = datetime(2021,05,01);
dry_end = datetime(2021,11,01);
dt = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Full Time Series

%%%%%%%%%%%%%%%%
%%% Salinity %%%
%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,~,s2x,~,s2y,k,~] = xcovar(oism.avgsalt_1dayrunmean, hatfield_HT.salt_interp, 10);

%%% Getting Pearson Correlation Coefficient from covariance
corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);

%%% Identifying max correlation coefficient
[M,I] = max(corr_coeff);
disp('   ');
disp('Highest Salinity Correlation Coefficient: ' + string(M));

%%% Calculating the corresponding time lag
disp('Corresponding Time Lag: ' + string(k(I)*dt) + ' Day(s)');

%%% Creating variables for parabolic fitting
idy = I-2:1:I+2;
corr_coeff_4_fitting = corr_coeff(idy);
time_lag_4_fitting = k(idy)*dt;

%%% Fitting parabola
p = polyfit(time_lag_4_fitting, corr_coeff_4_fitting, 2);

%%% Finding max
fitted_time_lag = roots(polyder(p));

disp('Fitted Time Lag: ' + string(fitted_time_lag) + ' Days(s)');

%%% Visualizing parabolic fitting
f = polyval(p, k(I)-2:0.1:k(I)+2);
figure()
plot(k, corr_coeff)
hold on
plot(k(I)-2:0.1:k(I)+2, f);
xline(fitted_time_lag)
hold off
xlabel('Time Lag')
ylabel('Correlation Coefficient')
title('Salinity')

clear Rxy mux s2x muy s2y k Nk M I corr_coeff corr_coeff_4_fitting idy time_lag_4_fitting fitted_time_lag f

%%%%%%%%%%%%%%%%%%%
%%% Temperature %%%
%%%%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,~,s2x,~,s2y,k,~] = xcovar(oism.avgtemp_1dayrunmean, hatfield_HT.temp_interp, 10);

%%% Getting Pearson Correlation Coefficient from covariance
corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);

%%% Identifying max correlation coefficient
[M,I] = max(corr_coeff);
disp('   ');
disp('Highest Temperature Correlation Coefficient: ' + string(M));

%%% Calculating the corresponding time lag
disp('Corresponding Time Lag: ' + string(k(I)*dt) + ' Day(s)');

%%% Creating variables for parabolic fitting
idy = I-2:1:I+2;
corr_coeff_4_fitting = corr_coeff(idy);
time_lag_4_fitting = k(idy)*dt;

%%% Fitting parabola
p = polyfit(time_lag_4_fitting, corr_coeff_4_fitting, 2);

%%% Finding max
fitted_time_lag = roots(polyder(p));

disp('Fitted Time Lag: ' + string(fitted_time_lag) + ' Days(s)');

%%% Visualizing parabolic fitting
f = polyval(p, k(I)-2:0.1:k(I)+2);
figure()
plot(k, corr_coeff)
hold on
plot(k(I)-2:0.1:k(I)+2, f);
xline(fitted_time_lag)
hold off
xlabel('Time Lag')
ylabel('Correlation Coefficient')
title('Temperature')

clear Rxy mux s2x muy s2y k Nk M I corr_coeff corr_coeff_4_fitting idy time_lag_4_fitting fitted_time_lag f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Dry Season

%%%%%%%%%%%%%%%%
%%% Salinity %%%
%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,~,s2x,~,s2y,k,~] = xcovar(oism.avgsalt_1dayrunmean(oism.datetime_1dayrunmean > dry_start & oism.datetime_1dayrunmean < dry_end), hatfield_HT.salt_interp(hatfield_HT.datetime_interp > dry_start & hatfield_HT.datetime_interp < dry_end), 10);

%%% Getting Pearson Correlation Coefficient from covariance
corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);

%%% Identifying max correlation coefficient
[M,I] = max(corr_coeff);
disp('   ');
disp('Highest Dry Season Salinity Correlation Coefficient: ' + string(M));

%%% Calculating the corresponding time lag
disp('Corresponding Time Lag: ' + string(k(I)*dt) + ' Day(s)');

%%% Creating variables for parabolic fitting
idy = I-2:1:I+2;
corr_coeff_4_fitting = corr_coeff(idy);
time_lag_4_fitting = k(idy)*dt;

%%% Fitting parabola
p = polyfit(time_lag_4_fitting, corr_coeff_4_fitting, 2);

%%% Finding max
fitted_time_lag = roots(polyder(p));

disp('Fitted Time Lag: ' + string(fitted_time_lag) + ' Days(s)');

%%% Visualizing parabolic fitting
f = polyval(p, k(I)-2:0.1:k(I)+2);
figure()
plot(k, corr_coeff)
hold on
plot(k(I)-2:0.1:k(I)+2, f);
xline(fitted_time_lag)
hold off
xlabel('Time Lag')
ylabel('Correlation Coefficient')
title('Salinity - Dry Season')

clear Rxy mux s2x muy s2y k Nk M I corr_coeff corr_coeff_4_fitting idy time_lag_4_fitting fitted_time_lag f

%%%%%%%%%%%%%%%%%%%
%%% Temperature %%%
%%%%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,~,s2x,~,s2y,k,~] = xcovar(oism.avgtemp_1dayrunmean(oism.datetime_1dayrunmean > dry_start & oism.datetime_1dayrunmean < dry_end), hatfield_HT.temp_interp(hatfield_HT.datetime_interp > dry_start & hatfield_HT.datetime_interp < dry_end), 10);

%%% Getting Pearson Correlation Coefficient from covariance
corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);

%%% Identifying max correlation coefficient
[M,I] = max(corr_coeff);
disp('   ');
disp('Highest Dry Season Temperature Correlation Coefficient: ' + string(M));

%%% Calculating the corresponding time lag
disp('Corresponding Time Lag: ' + string(k(I)*dt) + ' Day(s)');

%%% Creating variables for parabolic fitting
idy = I-2:1:I+2;
corr_coeff_4_fitting = corr_coeff(idy);
time_lag_4_fitting = k(idy)*dt;

%%% Fitting parabola
p = polyfit(time_lag_4_fitting, corr_coeff_4_fitting, 2);

%%% Finding max
fitted_time_lag = roots(polyder(p));

disp('Fitted Time Lag: ' + string(fitted_time_lag) + ' Days(s)');

%%% Visualizing parabolic fitting
f = polyval(p, k(I)-2:0.1:k(I)+2);
figure()
plot(k, corr_coeff)
hold on
plot(k(I)-2:0.1:k(I)+2, f);
xline(fitted_time_lag)
hold off
xlabel('Time Lag')
ylabel('Correlation Coefficient')
title('Temperature - Dry Season')

clear Rxy mux s2x muy s2y k Nk M I corr_coeff corr_coeff_4_fitting idy time_lag_4_fitting fitted_time_lag f

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Wet Season

%%%%%%%%%%%%%%%%
%%% Salinity %%%
%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,~,s2x,~,s2y,k,~] = xcovar(oism.avgsalt_1dayrunmean(oism.datetime_1dayrunmean < dry_start | oism.datetime_1dayrunmean > dry_end), hatfield_HT.salt_interp(hatfield_HT.datetime_interp < dry_start | hatfield_HT.datetime_interp > dry_end), 10);

%%% Getting Pearson Correlation Coefficient from covariance
corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);

%%% Identifying max correlation coefficient
[M,I] = max(corr_coeff);
disp('   ');
disp('Highest Wet Season Salinity Correlation Coefficient: ' + string(M));

%%% Calculating the corresponding time lag
disp('Corresponding Time Lag: ' + string(k(I)*dt) + ' Day(s)');

%%% Creating variables for parabolic fitting
idy = I-2:1:I+2;
corr_coeff_4_fitting = corr_coeff(idy);
time_lag_4_fitting = k(idy)*dt;

%%% Fitting parabola
p = polyfit(time_lag_4_fitting, corr_coeff_4_fitting, 2);

%%% Finding max
fitted_time_lag = roots(polyder(p));

disp('Fitted Time Lag: ' + string(fitted_time_lag) + ' Days(s)');

%%% Visualizing parabolic fitting
f = polyval(p, k(I)-2:0.1:k(I)+2);
figure()
plot(k, corr_coeff)
hold on
plot(k(I)-2:0.1:k(I)+2, f);
xline(fitted_time_lag)
hold off
xlabel('Time Lag')
ylabel('Correlation Coefficient')
title('Salinity - Wet Season')

clear Rxy mux s2x muy s2y k Nk M I corr_coeff corr_coeff_4_fitting idy time_lag_4_fitting fitted_time_lag f

%%%%%%%%%%%%%%%%%%%
%%% Temperature %%%
%%%%%%%%%%%%%%%%%%%

%%% Calculating covariance
[Rxy,~,s2x,~,s2y,k,~] = xcovar(oism.avgtemp_1dayrunmean(oism.datetime_1dayrunmean < dry_start | oism.datetime_1dayrunmean > dry_end), hatfield_HT.temp_interp(hatfield_HT.datetime_interp < dry_start | hatfield_HT.datetime_interp > dry_end), 10);

%%% Getting Pearson Correlation Coefficient from covariance
corr_coeff = Rxy./((s2x).^0.5.*(s2y).^0.5);

%%% Identifying max correlation coefficient
[M,I] = max(corr_coeff);
disp('   ');
disp('Highest Wet Season Temperature Correlation Coefficient: ' + string(M));

%%% Calculating the corresponding time lag
disp('Corresponding Time Lag: ' + string(k(I)*dt) + ' Day(s)');

%%% Creating variables for parabolic fitting
idy = I-2:1:I+2;
corr_coeff_4_fitting = corr_coeff(idy);
time_lag_4_fitting = k(idy)*dt;

%%% Fitting parabola
p = polyfit(time_lag_4_fitting, corr_coeff_4_fitting, 2);

%%% Finding max
fitted_time_lag = roots(polyder(p));

disp('Fitted Time Lag: ' + string(fitted_time_lag) + ' Days(s)');

%%% Visualizing parabolic fitting
f = polyval(p, k(I)-2:0.1:k(I)+2);
figure()
plot(k, corr_coeff)
hold on
plot(k(I)-2:0.1:k(I)+2, f);
xline(fitted_time_lag)
hold off
xlabel('Time Lag')
ylabel('Correlation Coefficient')
title('Temperature - Wet Season')

clear Rxy mux s2x muy s2y k Nk M I corr_coeff corr_coeff_4_fitting idy time_lag_4_fitting fitted_time_lag f
