%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Loading Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Clearing Workspace
clear

%%% Loading OISM Data
load('OOI_CE_OISM');

%%% Loading OSSM Data
load('OOI_CE_OSSM');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Computing Running Mean of Data %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finding the first and last days of the time series
first_day = datenum(min(oism.salt7m.datetime));
last_day = datenum(max(oism.salt7m.datetime));

%%% Creating an hourly time grid
timegrid = floor(first_day):datenum(hours(1)):ceil(last_day);

%%%%%%%%%%%%%%%%%
%%% OISM DATA %%%
%%%%%%%%%%%%%%%%%

oism.datenum_rm1hr = timegrid;
oism.datetime_rm1hr = datetime(timegrid, 'ConvertFrom', 'datenum');

%%% 7m runmean for salinity, temperature, pressure, and DO
oism.salt7m_rm1hr = runmean(timegrid, datenum(oism.salt7m.datetime), oism.salt7m.salt);
oism.temp7m_rm1hr = runmean(timegrid, datenum(oism.temp7m.datetime), oism.temp7m.temp);
oism.pres7m_rm1hr = runmean(timegrid, datenum(oism.pres7m.datetime), oism.pres7m.pres);
oism.do7m_rm1hr = runmean(timegrid, datenum(oism.do7m.datetime), oism.do7m.do);

%%% 25m runmean for salinity, temperature, pressure, and DO
oism.salt25m_rm1hr = runmean(timegrid, datenum(oism.salt25m.datetime), oism.salt25m.salt);
oism.temp25m_rm1hr = runmean(timegrid, datenum(oism.temp25m.datetime), oism.temp25m.temp);
oism.pres25m_rm1hr = runmean(timegrid, datenum(oism.pres25m.datetime), oism.pres25m.pres);
oism.do25m_rm1hr = runmean(timegrid, datenum(oism.do25m.datetime), oism.do25m.do);

%%%%%%%%%%%%%%%%%
%%% OSSM DATA %%%
%%%%%%%%%%%%%%%%%

ossm.datenum_rm1hr = timegrid;
ossm.datetime_rm1hr = datetime(timegrid, 'ConvertFrom', 'datenum');

%%% 7m runmean for salinity, temperature, and pressure
ossm.salt7m_rm1hr = runmean(timegrid, datenum(ossm.salt7m.datetime), ossm.salt7m.salt);
ossm.temp7m_rm1hr = runmean(timegrid, datenum(ossm.temp7m.datetime), ossm.temp7m.temp);
ossm.pres7m_rm1hr = runmean(timegrid, datenum(ossm.pres7m.datetime), ossm.pres7m.pres);
ossm.do7m_rm1hr = runmean(timegrid, datenum(ossm.do7m.datetime), ossm.do7m.do);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Computing % Saturation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%%% OISM DATA %%%
%%%%%%%%%%%%%%%%%

%%% Calculating Absolute Salinity
oism.abs_salt7m_rm1hr = gsw_SA_from_SP(oism.salt7m_rm1hr, oism.pres7m_rm1hr, oism.lon, oism.lat);
oism.abs_salt25m_rm1hr = gsw_SA_from_SP(oism.salt25m_rm1hr, oism.pres25m_rm1hr, oism.lon, oism.lat);

%%% Calculating potential temperature
oism.pot_temp7m_rm1hr = gsw_pt_from_t(oism.abs_salt7m_rm1hr, oism.temp7m_rm1hr, oism.pres7m_rm1hr);
oism.pot_temp25m_rm1hr = gsw_pt_from_t(oism.abs_salt25m_rm1hr, oism.temp25m_rm1hr, oism.pres25m_rm1hr);

%%% Calculating solubility of oxygen
oism.O2sol7m_rm1hr = gsw_O2sol_SP_pt(oism.salt7m_rm1hr, oism.pot_temp7m_rm1hr);
oism.O2sol25m_rm1hr = gsw_O2sol_SP_pt(oism.salt25m_rm1hr, oism.pot_temp25m_rm1hr);

%%% Calculating % Saturation
oism.pctO2sat7m_rm1hr = (oism.do7m_rm1hr ./ oism.O2sol7m_rm1hr) .* 100;
oism.pctO2sat25m_rm1hr = (oism.do25m_rm1hr ./ oism.O2sol25m_rm1hr) .* 100;

%%%%%%%%%%%%%%%%%
%%% OSSM DATA %%%
%%%%%%%%%%%%%%%%%

%%% Calculating Absolute Salinity
ossm.abs_salt7m_rm1hr = gsw_SA_from_SP(ossm.salt7m_rm1hr, ossm.pres7m_rm1hr, ossm.lon, ossm.lat);

%%% Calculating potential temperature
ossm.pot_temp7m_rm1hr = gsw_pt_from_t(ossm.abs_salt7m_rm1hr, ossm.temp7m_rm1hr, ossm.pres7m_rm1hr);

%%% Calculating solubility of oxygen
ossm.O2sol7m_rm1hr = gsw_O2sol_SP_pt(ossm.salt7m_rm1hr, ossm.pot_temp7m_rm1hr);

%%% Calculating % Saturation
ossm.pctO2sat7m_rm1hr = (ossm.do7m_rm1hr ./ ossm.O2sol7m_rm1hr) .* 100;

%%

figure('Renderer', 'painters', 'Position', [100 100 1200 800])

ax1 = subplot(311);
hold on
plot(oism.datetime_rm1hr, oism.do7m_rm1hr, 'k', 'DisplayName', 'Inshore - 7m')
plot(oism.datetime_rm1hr, oism.do25m_rm1hr, 'b', 'DisplayName', 'Inshore - 25m')
plot(ossm.datetime_rm1hr, ossm.do7m_rm1hr, 'r', 'DisplayName', 'Shelf - 7m')
hold off
ylabel('umol/kg')
title('Dissolved Oxygen Concentration');
legend()

ax2 = subplot(312);
hold on
plot(oism.datetime_rm1hr, oism.O2sol7m_rm1hr, 'k', 'DisplayName', 'Inshore - 7m')
plot(oism.datetime_rm1hr, oism.O2sol25m_rm1hr, 'b', 'DisplayName', 'Inshore - 25m')
plot(ossm.datetime_rm1hr, ossm.O2sol7m_rm1hr, 'r', 'DisplayName', 'Shelf - 7m')
hold off
ylabel('umol/kg')
title('Equilibrium Oxygen Concentration');
legend()

ax3 = subplot(313);
hold on
plot(oism.datetime_rm1hr, oism.pctO2sat7m_rm1hr, 'k', 'DisplayName', 'Inshore - 7m')
plot(oism.datetime_rm1hr, oism.pctO2sat25m_rm1hr, 'b', 'DisplayName', 'Inshore - 25m')
plot(ossm.datetime_rm1hr, ossm.pctO2sat7m_rm1hr, 'r', 'DisplayName', 'Shelf - 7m')
hold off
title('Percent Saturation');
legend()

linkaxes([ax1 ax2 ax3], 'x');


