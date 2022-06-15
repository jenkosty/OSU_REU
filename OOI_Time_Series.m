%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Loading and pre-processing OOI data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading OISM data at different depths
oism.temp1m = readtable('OOI_CE_OISM_SB_Temperature.csv');
oism.salt1m = readtable('OOI_CE_OISM_SB_Salinity.csv');
oism.temp7m = readtable('OOI_CE_OISM_NSIF_Temperature.csv');
oism.salt7m = readtable('OOI_CE_OISM_NSIF_Salinity.csv');
oism.temp25m = readtable('OOI_CE_OISM_SMFN_Temperature.csv');
oism.salt25m = readtable('OOI_CE_OISM_SMFN_Salinity.csv');

%%% Loading OOSM data at different depths
oosm.temp1m = readtable('OOI_CE_OOSM_SB_Temperature.csv');
oosm.salt1m = readtable('OOI_CE_OOSM_SB_Salinity.csv');
oosm.temp7m = readtable('OOI_CE_OOSM_NSIF_Temperature.csv');
oosm.salt7m = readtable('OOI_CE_OOSM_NSIF_Salinity.csv');

%%% Filtering OISM to the highest quality control indicators
oism.temp7m = oism.temp7m(oism.temp7m.sea_water_temperature_qc_agg == 1, :);
oism.salt7m = oism.salt7m(oism.salt7m.sea_water_practical_salinity_qc_agg == 1, :);
oism.temp25m = oism.temp25m(oism.temp25m.sea_water_temperature_qc_agg == 1, :);
oism.salt25m = oism.salt25m(oism.salt25m.sea_water_practical_salinity_qc_agg == 1, :);

%%% Converting time strings to datetime
oism.temp7m.datetime = datetime(vertcat(oism.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt7m.datetime = datetime(vertcat(oism.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.temp25m.datetime = datetime(vertcat(oism.temp25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt25m.datetime = datetime(vertcat(oism.salt25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Interpolating all data to a consistent time grid %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating a consistently spaced time grid
time_grid = (datetime(2014,04,01):hours(2):datetime(2022,06,01))';

%%% Interpolating onto time grid
oism_interp.time = time_grid;
oism_interp.temp7m = interp1gap(datenum(oism.temp7m.datetime), oism.temp7m.sea_water_temperature, datenum(time_grid), datenum(hours(3)));
oism_interp.salt7m = interp1gap(datenum(oism.salt7m.datetime), oism.salt7m.sea_water_practical_salinity, datenum(time_grid), datenum(hours(3)));
oism_interp.temp25m = interp1gap(datenum(oism.temp25m.datetime), oism.temp25m.sea_water_temperature, datenum(time_grid), datenum(hours(3)));
oism_interp.salt25m = interp1gap(datenum(oism.salt25m.datetime), oism.salt25m.sea_water_practical_salinity, datenum(time_grid), datenum(hours(3)));

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
%%% Plotting temperature and salinity time series (both raw and
%%% interpolated)

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



