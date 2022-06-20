=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating a consistently spaced time grid
time_grid = (datetime(2014,04,01):hours(2):datetime(2022,06,01))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Loading and pre-processing OOI data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OISM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oism.temp1m = readtable('OOI_CE_OISM_SB_Temperature.csv');
oism.salt1m = readtable('OOI_CE_OISM_SB_Salinity.csv');
oism.temp7m = readtable('OOI_CE_OISM_NSIF_Temperature.csv');
oism.salt7m = readtable('OOI_CE_OISM_NSIF_Salinity.csv');
oism.temp25m = readtable('OOI_CE_OISM_SMFN_Temperature.csv');
oism.salt25m = readtable('OOI_CE_OISM_SMFN_Salinity.csv');

%%% Filtering OISM to the highest quality control indicators
oism.temp1m = oism.temp1m(oism.temp1m.sea_water_temperature_qc_agg == 1, :);
oism.salt1m = oism.salt1m(oism.salt1m.sea_water_practical_salinity_qc_agg == 1, :);
oism.temp7m = oism.temp7m(oism.temp7m.sea_water_temperature_qc_agg == 1, :);
oism.salt7m = oism.salt7m(oism.salt7m.sea_water_practical_salinity_qc_agg == 1, :);
oism.temp25m = oism.temp25m(oism.temp25m.sea_water_temperature_qc_agg == 1, :);
oism.salt25m = oism.salt25m(oism.salt25m.sea_water_practical_salinity_qc_agg == 1, :);

%%% Converting time strings to datetime
oism.temp1m.datetime = datetime(vertcat(oism.temp1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt1m.datetime = datetime(vertcat(oism.salt1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.temp7m.datetime = datetime(vertcat(oism.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt7m.datetime = datetime(vertcat(oism.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.temp25m.datetime = datetime(vertcat(oism.temp25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt25m.datetime = datetime(vertcat(oism.salt25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
oism.temp1m = renamevars(oism.temp1m,'sea_water_temperature','temp');
oism.salt1m = renamevars(oism.salt1m,'sea_water_practical_salinity','salt');
oism.temp7m = renamevars(oism.temp7m,'sea_water_temperature','temp');
oism.salt7m = renamevars(oism.salt7m,'sea_water_practical_salinity','salt');
oism.temp25m = renamevars(oism.temp25m,'sea_water_temperature','temp');
oism.salt25m = renamevars(oism.salt25m,'sea_water_practical_salinity','salt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OSSM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ossm.temp1m = readtable('OOI_CE_OOSM_SB_Temperature.csv');
ossm.salt1m = readtable('OOI_CE_OOSM_SB_Salinity.csv');
ossm.temp7m = readtable('OOI_CE_OOSM_NSIF_Temperature.csv');
ossm.salt7m = readtable('OOI_CE_OOSM_NSIF_Salinity.csv');

%%% Filtering OSSM to the highest quality control indicators
ossm.temp1m = ossm.temp1m(ossm.temp1m.sea_surface_temperature_qc_agg == 1, :);
ossm.salt1m = ossm.salt1m(ossm.salt1m.sea_water_practical_salinity_qc_agg == 1, :);
ossm.temp7m = ossm.temp7m(ossm.temp7m.sea_water_temperature_qc_agg == 1, :);
ossm.salt7m = ossm.salt7m(ossm.salt7m.sea_water_practical_salinity_qc_agg == 1, :);

%%% Converting time strings to datetime
ossm.temp1m.datetime = datetime(vertcat(ossm.temp1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.salt1m.datetime = datetime(vertcat(ossm.salt1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.temp7m.datetime = datetime(vertcat(ossm.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.salt7m.datetime = datetime(vertcat(ossm.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
ossm.temp1m = renamevars(ossm.temp1m,'sea_surface_temperature','temp');
ossm.salt1m = renamevars(ossm.salt1m,'sea_water_practical_salinity','salt');
ossm.temp7m = renamevars(ossm.temp7m,'sea_water_temperature','temp');
ossm.salt7m = renamevars(ossm.salt7m,'sea_water_practical_salinity','salt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OOSM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oosm.temp1m = readtable('OOI_CE_OOSM_SB_Temperature.csv');
oosm.salt1m = readtable('OOI_CE_OOSM_SB_Salinity.csv');
oosm.temp7m = readtable('OOI_CE_OOSM_NSIF_Temperature.csv');
oosm.salt7m = readtable('OOI_CE_OOSM_NSIF_Salinity.csv');

%%% Filtering OOSM to the highest quality control indicators
oosm.temp1m = oosm.temp1m(oosm.temp1m.sea_surface_temperature_qc_agg == 1, :);
oosm.salt1m = oosm.salt1m(oosm.salt1m.sea_water_practical_salinity_qc_agg == 1, :);
oosm.temp7m = oosm.temp7m(oosm.temp7m.sea_water_temperature_qc_agg == 1, :);
oosm.salt7m = oosm.salt7m(oosm.salt7m.sea_water_practical_salinity_qc_agg == 1, :);

%%% Converting time strings to datetime
oosm.temp1m.datetime = datetime(vertcat(oosm.temp1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.salt1m.datetime = datetime(vertcat(oosm.salt1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.temp7m.datetime = datetime(vertcat(oosm.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.salt7m.datetime = datetime(vertcat(oosm.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
oosm.temp1m = renamevars(oosm.temp1m,'sea_surface_temperature','temp');
oosm.salt1m = renamevars(oosm.salt1m,'sea_water_practical_salinity','salt');
oosm.temp7m = renamevars(oosm.temp7m,'sea_water_temperature','temp');
oosm.salt7m = renamevars(oosm.salt7m,'sea_water_practical_salinity','salt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Loading and pre-processing meterological data %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading meterological data from different years
metero_unfmt.data2016 = readtable('2016meterologicaldata.txt');
metero_unfmt.data2017 = readtable('2017meterologicaldata.txt');
metero_unfmt.data2018 = readtable('2018meterologicaldata.txt');
metero_unfmt.data2019 = readtable('2019meterologicaldata.txt');
metero_unfmt.data2020 = readtable('2020meterologicaldata.txt');
metero_unfmt.data2021 = readtable('2021meterologicaldata.txt');

%%% Extracting wind direction data
metero.wind_dir = table2array(vertcat(metero_unfmt.data2016(:,6),metero_unfmt.data2017(:,6),metero_unfmt.data2018(:,6),...
    metero_unfmt.data2019(:,6), metero_unfmt.data2020(:,6), metero_unfmt.data2021(:,6)));

%%% Extracting wind speed data
metero.wind_spd = table2array(vertcat(metero_unfmt.data2016(:,7),metero_unfmt.data2017(:,7),metero_unfmt.data2018(:,7),...
    metero_unfmt.data2019(:,7), metero_unfmt.data2020(:,7), metero_unfmt.data2021(:,7)));

%%% Extracting time of measurement data
yrs = table2array(vertcat(metero_unfmt.data2016(:,1),metero_unfmt.data2017(:,1),metero_unfmt.data2018(:,1),...
    metero_unfmt.data2019(:,1), metero_unfmt.data2020(:,1), metero_unfmt.data2021(:,1)));

mths = table2array(vertcat(metero_unfmt.data2016(:,2),metero_unfmt.data2017(:,2),metero_unfmt.data2018(:,2),...
    metero_unfmt.data2019(:,2), metero_unfmt.data2020(:,2), metero_unfmt.data2021(:,2)));

days = table2array(vertcat(metero_unfmt.data2016(:,3),metero_unfmt.data2017(:,3),metero_unfmt.data2018(:,3),...
    metero_unfmt.data2019(:,3), metero_unfmt.data2020(:,3), metero_unfmt.data2021(:,3)));

hrs = table2array(vertcat(metero_unfmt.data2016(:,4),metero_unfmt.data2017(:,4),metero_unfmt.data2018(:,4),...
    metero_unfmt.data2019(:,4), metero_unfmt.data2020(:,4), metero_unfmt.data2021(:,4)));

mnts = table2array(vertcat(metero_unfmt.data2016(:,5),metero_unfmt.data2017(:,5),metero_unfmt.data2018(:,5),...
    metero_unfmt.data2019(:,5), metero_unfmt.data2020(:,5), metero_unfmt.data2021(:,5)));

snds = zeros(size(yrs)); %%% Setting seconds to 0

%%% Converting time of measurement data to datetime
metero.datetime = datetime(yrs,mths,days,hrs,mnts, snds);

clear yrs mths days hrs mnts snds metero_unfmt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Loading and pre-processing river discharge data %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading river discharge data
riverflow_unfmt = readtable('Station_14306030_instantaneous_flow.txt','Delimiter','\t');

%%% Extracting river discharge data. Converting from cubic feet to cubic
%%% meters per second. Scaling by 2.53 to account for additional inflow 
%%% downstream of gauge and upstream of estuary
riverflow.flow = riverflow_unfmt{:,3}*2.52*(.3048)^3;

%%% Extrating time of measurement data
riverflow.datetime = datetime(riverflow_unfmt{:,2}, 'InputFormat', 'MM-dd-yyyy HH:mm');

clear riverflow_unfmt

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Interpolating all data to a consistent time grid %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Setting the time max gap between adjacent measurements for allowed
%%% interpolation
max_gap = hours(3);

%%% Interpolating OISM data onto time grid
oism.time_interp = time_grid;
oism.temp1m_interp = interp1gap(datenum(oism.temp1m.datetime), oism.temp1m.temp, datenum(time_grid), datenum(max_gap));
oism.salt1m_interp = interp1gap(datenum(oism.salt1m.datetime), oism.salt1m.salt, datenum(time_grid), datenum(max_gap));
oism.temp7m_interp = interp1gap(datenum(oism.temp7m.datetime), oism.temp7m.temp, datenum(time_grid), datenum(max_gap));
oism.salt7m_interp = interp1gap(datenum(oism.salt7m.datetime), oism.salt7m.salt, datenum(time_grid), datenum(max_gap));
oism.temp25m_interp = interp1gap(datenum(oism.temp25m.datetime), oism.temp25m.temp, datenum(time_grid), datenum(max_gap));
oism.salt25m_interp = interp1gap(datenum(oism.salt25m.datetime), oism.salt25m.salt, datenum(time_grid), datenum(max_gap));

%%% Interpolating OOSM data onto time grid
ossm.time_interp = time_grid;
ossm.temp1m_interp = interp1gap(datenum(ossm.temp1m.datetime), ossm.temp1m.temp, datenum(time_grid), datenum(max_gap));
ossm.salt1m_interp = interp1gap(datenum(ossm.salt1m.datetime), ossm.salt1m.salt, datenum(time_grid), datenum(max_gap));
ossm.temp7m_interp = interp1gap(datenum(ossm.temp7m.datetime), ossm.temp7m.temp, datenum(time_grid), datenum(max_gap));
ossm.salt7m_interp = interp1gap(datenum(ossm.salt7m.datetime), ossm.salt7m.salt, datenum(time_grid), datenum(max_gap));

%%% Interpolating OOSM data onto time grid
oosm.time_interp = time_grid;
oosm.temp1m_interp = interp1gap(datenum(oosm.temp1m.datetime), oosm.temp1m.temp, datenum(time_grid), datenum(max_gap));
oosm.salt1m_interp = interp1gap(datenum(oosm.salt1m.datetime), oosm.salt1m.salt, datenum(time_grid), datenum(max_gap));
oosm.temp7m_interp = interp1gap(datenum(oosm.temp7m.datetime), oosm.temp7m.temp, datenum(time_grid), datenum(max_gap));
oosm.salt7m_interp = interp1gap(datenum(oosm.salt7m.datetime), oosm.salt7m.salt, datenum(time_grid), datenum(max_gap));

%%% Interpolating wind data onto time grid
metero.time_interp = time_grid;
metero.wind_dir_interp = interp1gap(datenum(metero.datetime), metero.wind_dir, datenum(time_grid), datenum(max_gap));
metero.wind_spd_interp = interp1gap(datenum(metero.datetime), metero.wind_spd, datenum(time_grid), datenum(max_gap));

%%% Interpolating river discharge data
riverflow.time_interp = time_grid;
riverflow.flow_interp = interp1gap(datenum(riverflow.datetime), riverflow.flow, datenum(time_grid), datenum(max_gap));
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Creating a running mean for all data %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Setting the time step for the running mean
runmean_time_step = datenum(hours(1));

%%% Creating running mean for the OISM data on time grid
oism.time_runmean = time_grid;
oism.temp1m_runmean = runmean(time_grid, runmean_time_step, oism.temp1m.datetime, oism.temp1m);
oism.salt1m_runmean = runmean(time_grid, runmean_time_step, oism.salt1m.datetime, oism.salt1m.sea_water_practical_salinity);
oism.temp7m_runmean = runmean(time_grid, runmean_time_step, oism.temp7m.datetime, oism.temp7m.temp);
oism.salt7m_runmean = runmean(time_grid, runmean_time_step, oism.salt7m.datetime, oism.salt7m.salt);
oism.temp25m_runmean = runmean(time_grid, runmean_time_step, oism.temp25m.datetime, oism.temp25m);
oism.salt25m_runmean = runmean(time_grid, runmean_time_step, oism.salt25m.datetime, oism.salt25m);

%%% Creating running mean for the OSSM data on time grid
ossm.time_runmean = time_grid;
ossm.temp1m_runmean = runmean(time_grid, runmean_time_step, ossm.temp1m.datetime, ossm.temp1m);
ossm.salt1m_runmean = runmean(time_grid, runmean_time_step, ossm.salt1m.datetime, ossm.salt1m);
ossm.temp7m_runmean = runmean(time_grid, runmean_time_step, ossm.temp7m.datetime, ossm.temp7m);
ossm.salt7m_runmean = runmean(time_grid, runmean_time_step, ossm.salt7m.datetime, ossm.salt7m);

%%% Creating running mean for the OOSM data on time grid
oosm.time_runmean = time_grid;
oosm.temp1m_runmean = runmean(time_grid, runmean_time_step, oosm.temp1m.datetime, oosm.temp1m);
oosm.salt1m_runmean = runmean(time_grid, runmean_time_step, oosm.salt1m.datetime, oosm.salt1m);
oosm.temp7m_runmean = runmean(time_grid, runmean_time_step, oosm.temp7m.datetime, oosm.temp7m);
oosm.salt7m_runmean = runmean(time_grid, runmean_time_step, oosm.salt7m.datetime, oosm.salt7m);

%%% Creating running mean for wind data on time grid
metero.time_runmean = time_grid;
metero.wind_dir_runmean = runmean(time_grid, runmean_time_step, metero.datetime, metero.wind_dir);
metero.wind_spd_runmean = runmean(time_grid, runmean_time_step, metero.datetime, metero.wind_spd);

%%% Creating running mean for river discharge data on time grid
riverflow.time_runmean = time_grid;
riverflow.flow_runmean = runmean(time_grid, runmean_time_step, riverflow.datetime, riverflow.flow);

%%
% %%% Plotting time gaps between measurements (used to determine appropriate
% %%% time grid for interpolation)
% 
% %%% OISM Near Surface Instrument Frame
% figure('Renderer', 'painters', 'Position', [100 100 1000 600])
% sgtitle('Oregon Inshoor Surface Mooring - Near Surface Instrument Frame (7 m) - Temperature')
% 
% %%% Full scatter plot
% subplot(211);
% scatter(oism.temp7m.datetime(1:end-1), hours(diff(vertcat(oism.temp7m.datetime))), 'filled');
% ylabel('Hours Between Measurements')
% 
% %%% Zooming in
% subplot(212);
% scatter(oism.temp7m.datetime(1:end-1), hours(diff(vertcat(oism.temp7m.datetime))), 'filled');
% ylabel('Hours Between Measurements')
% ylim([0 5])
% 
% %%% OSIM Seafloor Multi-Function Node
% figure('Renderer', 'painters', 'Position', [100 100 1000 600])
% sgtitle('Oregon Inshoor Surface Mooring - Seafloor Multi-Function Node (25 m) - Temperature')
% 
% %%% Full scatter plot
% subplot(211);
% scatter(oism.temp25m.datetime(1:end-1), hours(diff(vertcat(oism.temp25m.datetime))), 'filled');
% ylabel('Hours Between Measurements')
% 
% %%% Zooming in
% subplot(212);
% scatter(oism.temp25m.datetime(1:end-1), hours(diff(vertcat(oism.temp25m.datetime))), 'filled');
% ylabel('Hours Between Measurements')
% ylim([0 5])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Plotting temperature and salinity time series (both raw and %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% interpolated) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% OISM Near Surface Instrument Frame
figure('Renderer', 'painters', 'Position', [100 100 1000 800])
sgtitle('Oregon Inshoor Surface Mooring - Near Surface Instrument Frame (7 m)');

    subplot(411)
    plot(oism_interp.time, oism_interp.temp7m) % interpolated temperature
    title('Interpolated')
    ylabel('Temperature (degC)');

    subplot(412)
    plot(oism.temp7m.datetime, oism.temp7m.sea_water_temperature) % raw temperature
    title('Raw')
    ylabel('Temperature (degC)');

    subplot(413)
    plot(oism_interp.time, oism_interp.salt7m) % interpolated salinity
    title('Interpolated')
    ylabel('Salinity (g/kg)');

    subplot(414)
    plot(oism.salt7m.datetime, oism.salt7m.sea_water_practical_salinity) % raw salinity
    title('Raw')
    ylabel('Salinity (g/kg)');

%%% OISM Seafloor Multi-Function Node
figure('Renderer', 'painters', 'Position', [100 100 1000 800])
sgtitle('Oregon Inshoor Surface Mooring - Seafloor Multi-Function Node (25 m)');

    subplot(411)
    plot(oism_interp.time, oism_interp.temp25m) % interpolated temperature
    title('Interpolated')
    ylabel('Temperature (degC)');

    subplot(412)
    plot(oism.temp25m.datetime, oism.temp25m.sea_water_temperature) % raw temperature
    title('Raw')
    ylabel('Temperature (degC)');

    subplot(413)
    plot(oism_interp.time, oism_interp.salt25m) % interpolated salinity
    title('Interpolated')
    ylabel('Salinity (g/kg)');

    subplot(414)
    plot(oism.salt25m.datetime, oism.salt25m.sea_water_practical_salinity) % raw salinity
    title('Raw')
    ylabel('Salinity (g/kg)');
    
%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Comparing different time series %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oism_var = oism_interp.salt25m;
oism_label = 'OISM 25m Salinity (g/kg)';

figure('Renderer', 'painters', 'Position', [100 100 1000 800])
sgtitle('OISM Salinity at 25m')

%%% OISM Variable vs OOSM 1m Salinity
subplot(311)
yyaxis left
plot(oism_interp.time, oism_var)
ylabel(oism_label);

yyaxis right
plot(oism_interp.time, oosm_interp.salt1m)
ylabel('OOSM 1m Salinity (g/kg)');

%%% OISM Variable vs Yaquina River discharge
subplot(312)
yyaxis left
plot(oism_interp.time, oism_var)
ylabel(oism_label)

yyaxis right
plot(riverflow_interp.time, riverflow_interp.flow)
ylabel('Yaquina River Discharge (m^3/s)');

%%% OISM Variable vs Wind Speed
subplot(313)
yyaxis left
plot(oism_interp.time, oism_var)
ylabel(oism_label)

yyaxis right
plot(metero_interp.time, metero_interp.wind_spd)
ylabel('Wind Speed (m/s)');
ylim([0 20])

clear oism_var oism_label
