%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
oism.pco2_7m = readtable('OOI_CE_OISM_NSIF_pCO2.csv');
oism.do7m = readtable('OOI_CE_OISM_NSIF_DO.csv');
oism.temp25m = readtable('OOI_CE_OISM_SMFN_Temperature.csv');
oism.salt25m = readtable('OOI_CE_OISM_SMFN_Salinity.csv');
oism.pco2_25m = readtable('OOI_CE_OISM_SMFN_pCO2.csv');
oism.do25m = readtable('OOI_CE_OISM_SMFN_DO.csv');

%%% Filtering OISM to the highest quality control indicators
%%% Note: DO and pCO2 data do not have quality control indicators
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
oism.pco2_7m.datetime = datetime(vertcat(oism.pco2_7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.do7m.datetime = datetime(vertcat(oism.do7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.temp25m.datetime = datetime(vertcat(oism.temp25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt25m.datetime = datetime(vertcat(oism.salt25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.pco2_25m.datetime = datetime(vertcat(oism.pco2_25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.do25m.datetime = datetime(vertcat(oism.do25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
oism.temp1m = renamevars(oism.temp1m,'sea_water_temperature','temp');
oism.salt1m = renamevars(oism.salt1m,'sea_water_practical_salinity','salt');
oism.temp7m = renamevars(oism.temp7m,'sea_water_temperature','temp');
oism.salt7m = renamevars(oism.salt7m,'sea_water_practical_salinity','salt');
oism.pco2_7m = renamevars(oism.pco2_7m,'partial_pressure_of_carbon_dioxide_in_sea_water','pco2');
oism.do7m = renamevars(oism.do7m,'moles_of_oxygen_per_unit_mass_in_sea_water','do');
oism.temp25m = renamevars(oism.temp25m,'sea_water_temperature','temp');
oism.salt25m = renamevars(oism.salt25m,'sea_water_practical_salinity','salt');
oism.pco2_25m = renamevars(oism.pco2_25m,'partial_pressure_of_carbon_dioxide_in_sea_water','pco2');
oism.do25m = renamevars(oism.do25m,'moles_of_oxygen_per_unit_mass_in_sea_water','do');

%%% Filtering pCO2 values to remove obviously bad data
oism.pco2_7m.pco2(oism.pco2_7m.pco2 > 1500) = NaN;
oism.pco2_25m.pco2(oism.pco2_25m.pco2 > 1500) = NaN;

%%% Removing bad pCO2 data (according to OOI annotations)
bad_data_7m = readtable('Bad_Data_OOI_CE_OISM_NSIF_pCO2.csv');
for i = 1:height(bad_data_7m)
    oism.pco2_7m.pco2(oism.pco2_7m.datetime > bad_data_7m.StartDate(i) & oism.pco2_7m.datetime < bad_data_7m.EndDate(i)) = NaN;
end

bad_data_25m = readtable('Bad_Data_OOI_CE_OISM_SMFN_pCO2.csv');
for i = 1:height(bad_data_25m)
    oism.pco2_25m.pco2(oism.pco2_25m.datetime > bad_data_25m.StartDate(i) & oism.pco2_25m.datetime < bad_data_25m.EndDate(i)) = NaN;
end

clear bad_data_7m bad_data_25m i

%%% Filtering DO values to remove obviously bad data
oism.do7m.do(oism.do7m.do <= 0) = NaN;
oism.do25m.do(oism.do25m.do <= 0) = NaN;

%%% Removing bad DO data (according to OOI annotations)
bad_data_7m = readtable('Bad_Data_OOI_CE_OISM_NSIF_DO.csv');
for i = 1:height(bad_data_7m)
    oism.do7m.do(oism.do7m.datetime > bad_data_7m.StartDate(i) & oism.do7m.datetime < bad_data_7m.EndDate(i)) = NaN;
end

bad_data_25m = readtable('Bad_Data_OOI_CE_OISM_SMFN_DO.csv');
for i = 1:height(bad_data_25m)
    oism.do25m.do(oism.do25m.datetime > bad_data_25m.StartDate(i) & oism.do25m.datetime < bad_data_25m.EndDate(i)) = NaN;
end

clear bad_data_7m bad_data_25m i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OISM QC'd DO Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading 7m QC DO Data 
load('OOI_CE_OISM_NSIF_QC_DO.mat')

%%% Formatting
do7m_qc.do = dissolved_oxygen';
do7m_qc.datetime = datetime(time,'ConvertFrom','datenum')';
do7m_qc.fail_biofoul_datetime = datetime(fail_biofoul, 'ConvertFrom', 'datenum')';

%%% Saving to structure
oism.do7m_qc = do7m_qc;

clear dissolved_oxygen fail_biofoul time do7m_qc

figure()
subplot(311)
plot(oism.do7m_qc.datetime, oism.do7m_qc.do)

oism.do7m_qc.do(ismember(oism.do7m_qc.datetime, oism.do7m_qc.fail_biofoul_datetime)) = NaN;

subplot(312)
plot(oism.do7m_qc.datetime, oism.do7m_qc.do)

subplot(313)
plot(oism.do7m.datetime, oism.do7m.do)

%%% Loading 25m QC DO Data
load('OOI_CE_OISM_SMFN_QC_DO.mat')

%%% Formatting
do25m_qc.do = dissolved_oxygen';
do25m_qc.datetime = datetime(time,'ConvertFrom','datenum')';
do25m_qc.fail_biofoul_datetime = datetime(fail_biofoul, 'ConvertFrom', 'datenum')';

%%% Saving to structure
oism.do25m_qc = do25m_qc;

clear dissolved_oxygen fail_biofoul time do25m_qc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OSSM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ossm.temp1m = readtable('OOI_CE_OSSM_SB_Temperature.csv');
ossm.salt1m = readtable('OOI_CE_OSSM_SB_Salinity.csv');
ossm.pco2_1m = readtable('OOI_CE_OSSM_SB_pCO2.csv');
ossm.temp7m = readtable('OOI_CE_OSSM_NSIF_Temperature.csv');
ossm.salt7m = readtable('OOI_CE_OSSM_NSIF_Salinity.csv');
ossm.do7m = readtable('OOI_CE_OSSM_NSIF_DO.csv');

%%% Filtering OSSM to the highest quality control indicators
%%% Note: DO and pCO2 data do not have quality control indicators
ossm.temp1m = ossm.temp1m(ossm.temp1m.sea_surface_temperature_qc_agg == 1, :);
ossm.salt1m = ossm.salt1m(ossm.salt1m.sea_water_practical_salinity_qc_agg == 1, :);
ossm.temp7m = ossm.temp7m(ossm.temp7m.sea_water_temperature_qc_agg == 1, :);
ossm.salt7m = ossm.salt7m(ossm.salt7m.sea_water_practical_salinity_qc_agg == 1, :);

%%% Converting time strings to datetime
ossm.temp1m.datetime = datetime(vertcat(ossm.temp1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.salt1m.datetime = datetime(vertcat(ossm.salt1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.pco2_1m.datetime = datetime(vertcat(ossm.pco2_1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.temp7m.datetime = datetime(vertcat(ossm.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.salt7m.datetime = datetime(vertcat(ossm.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.do7m.datetime = datetime(vertcat(ossm.do7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
ossm.temp1m = renamevars(ossm.temp1m,'sea_surface_temperature','temp');
ossm.salt1m = renamevars(ossm.salt1m,'sea_water_practical_salinity','salt');
ossm.pco2_1m = renamevars(ossm.pco2_1m,'surface_partial_pressure_of_carbon_dioxide_in_sea_water','pco2');
ossm.temp7m = renamevars(ossm.temp7m,'sea_water_temperature','temp');
ossm.salt7m = renamevars(ossm.salt7m,'sea_water_practical_salinity','salt');
ossm.do7m = renamevars(ossm.do7m,'moles_of_oxygen_per_unit_mass_in_sea_water','do');

%%% Filtering pCO2 values to remove obviously bad data
ossm.pco2_1m.pco2(ossm.pco2_1m.pco2 > 1500) = NaN;

%%% Removing bad pCO2 data (according to OOI annotations)
bad_data_1m = readtable('Bad_Data_OOI_CE_OSSM_SB_pCO2.csv');
for i = 1:height(bad_data_1m)
    ossm.pco2_1m.pco2(ossm.pco2_1m.datetime > bad_data_1m.StartDate(i) & ossm.pco2_1m.datetime < bad_data_1m.EndDate(i)) = NaN;
end

%%% Filtering DO values to remove obviously bad data
ossm.do7m.do(ossm.do7m.do < 0) = NaN;

%%% Removing bad DO data (according to OOI annotations)
bad_data_7m = readtable('Bad_Data_OOI_CE_OSSM_NSIF_DO.csv');
for i = 1:height(bad_data_7m)
    ossm.do7m.do(ossm.do7m.datetime > bad_data_7m.StartDate(i) & ossm.do7m.datetime < bad_data_7m.EndDate(i)) = NaN;
end

clear bad_data_1m bad_data_7m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OOSM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oosm.temp1m = readtable('OOI_CE_OOSM_SB_Temperature.csv');
oosm.salt1m = readtable('OOI_CE_OOSM_SB_Salinity.csv');
oosm.pco2_1m = readtable('OOI_CE_OOSM_SB_pCO2.csv');
oosm.temp7m = readtable('OOI_CE_OOSM_NSIF_Temperature.csv');
oosm.salt7m = readtable('OOI_CE_OOSM_NSIF_Salinity.csv');
oosm.do7m = readtable('OOI_CE_OOSM_NSIF_DO.csv');

%%% Filtering OOSM to the highest quality control indicators
%%% Note: DO and pCO2 data do not have quality control indicators
oosm.temp1m = oosm.temp1m(oosm.temp1m.sea_surface_temperature_qc_agg == 1, :);
oosm.salt1m = oosm.salt1m(oosm.salt1m.sea_water_practical_salinity_qc_agg == 1, :);
oosm.temp7m = oosm.temp7m(oosm.temp7m.sea_water_temperature_qc_agg == 1, :);
oosm.salt7m = oosm.salt7m(oosm.salt7m.sea_water_practical_salinity_qc_agg == 1, :);

%%% Converting time strings to datetime
oosm.temp1m.datetime = datetime(vertcat(oosm.temp1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.salt1m.datetime = datetime(vertcat(oosm.salt1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.pco2_1m.datetime = datetime(vertcat(oosm.pco2_1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.temp7m.datetime = datetime(vertcat(oosm.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.salt7m.datetime = datetime(vertcat(oosm.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oosm.do7m.datetime = datetime(vertcat(oosm.do7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
oosm.temp1m = renamevars(oosm.temp1m,'sea_surface_temperature','temp');
oosm.salt1m = renamevars(oosm.salt1m,'sea_water_practical_salinity','salt');
oosm.pco2_1m = renamevars(oosm.pco2_1m,'surface_partial_pressure_of_carbon_dioxide_in_sea_water','pco2');
oosm.temp7m = renamevars(oosm.temp7m,'sea_water_temperature','temp');
oosm.salt7m = renamevars(oosm.salt7m,'sea_water_practical_salinity','salt');
oosm.do7m = renamevars(oosm.do7m,'moles_of_oxygen_per_unit_mass_in_sea_water','do');

%%% Filtering pCO2 values to remove obviously bad data
oosm.pco2_1m.pco2(oosm.pco2_1m.pco2 > 1500) = NaN;

%%% Removing bad pCO2 data (according to OOI annotations)

%%% Filtering DO values to remove obviously bad data
oosm.do7m.do(oosm.do7m.do < 0) = NaN;

%%% Removing bad DO data (according to OOI annotations)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Loading and pre-processing meterological data %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%
%%% Shelf Data %%%
%%%%%%%%%%%%%%%%%%

%%% Loading meterological data from different years
metero_unfmt.data2016 = readtable('Shelf Wind Data/2016meterologicaldata.txt');
metero_unfmt.data2017 = readtable('Shelf Wind Data/2017meterologicaldata.txt');
metero_unfmt.data2018 = readtable('Shelf Wind Data/2018meterologicaldata.txt');
metero_unfmt.data2019 = readtable('Shelf Wind Data/2019meterologicaldata.txt');
metero_unfmt.data2020 = readtable('Shelf Wind Data/2020meterologicaldata.txt');
metero_unfmt.data2021 = readtable('Shelf Wind Data/2021meterologicaldata.txt');

%%% Extracting wind direction data
metero_shelf.wind_dir = table2array(vertcat(metero_unfmt.data2016(:,6),metero_unfmt.data2017(:,6),metero_unfmt.data2018(:,6),...
    metero_unfmt.data2019(:,6), metero_unfmt.data2020(:,6), metero_unfmt.data2021(:,6)));

%%% Extracting wind speed data
metero_shelf.wind_spd = table2array(vertcat(metero_unfmt.data2016(:,7),metero_unfmt.data2017(:,7),metero_unfmt.data2018(:,7),...
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
metero_shelf.datetime = datetime(yrs,mths,days,hrs,mnts,snds);

%%% Removing missing data
bad_data = find(metero_shelf.wind_spd > 90 | metero_shelf.wind_dir > 900);
metero_shelf.wind_dir(bad_data) = [];
metero_shelf.wind_spd(bad_data) = [];
metero_shelf.datetime(bad_data) = [];

clear yrs mths days hrs mnts snds metero_unfmt bad_data

%%%%%%%%%%%%%%%%%%%%%
%%% Offshore Data %%%
%%%%%%%%%%%%%%%%%%%%%

%%% Loading meterological data from different years
metero_unfmt.data2016 = readtable('Offshore Wind Data/2016meterologicaldata.txt');
metero_unfmt.data2017 = readtable('Offshore Wind Data/2017meterologicaldata.txt');
metero_unfmt.data2018 = readtable('Offshore Wind Data/2018meterologicaldata.txt');
metero_unfmt.data2019 = readtable('Offshore Wind Data/2019meterologicaldata.txt');
metero_unfmt.data2020 = readtable('Offshore Wind Data/2020meterologicaldata.txt');
metero_unfmt.data2021 = readtable('Offshore Wind Data/2021meterologicaldata.txt');

%%% Extracting wind direction data
metero_offshore.wind_dir = table2array(vertcat(metero_unfmt.data2016(:,6),metero_unfmt.data2017(:,6),metero_unfmt.data2018(:,6),...
    metero_unfmt.data2019(:,6), metero_unfmt.data2020(:,6), metero_unfmt.data2021(:,6)));

%%% Extracting wind speed data
metero_offshore.wind_spd = table2array(vertcat(metero_unfmt.data2016(:,7),metero_unfmt.data2017(:,7),metero_unfmt.data2018(:,7),...
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
metero_offshore.datetime = datetime(yrs,mths,days,hrs,mnts, snds);

%%% Removing missing data
bad_data = find(metero_offshore.wind_spd > 90 | metero_offshore.wind_dir > 900);
metero_offshore.wind_dir(bad_data) = [];
metero_offshore.wind_spd(bad_data) = [];
metero_offshore.datetime(bad_data) = [];

clear yrs mths days hrs mnts snds metero_unfmt bad_data

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
%%%%%%%%%%%% Loading and pre-processing Yaquina station data %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading station data
load('HFhightide.mat')

%%% Formatting in a structure
yaquina_station.datetime = tm;
yaquina_station.temp = T;
yaquina_station.temp_std = Tsd;
yaquina_station.salt = S;
yaquina_station.salt_std = Ssd;

clear S Ssd T tm Tsd 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Interpolating all data to a consistent time grid %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Setting the time max gap between adjacent measurements for allowed
%%% interpolation
max_gap = hours(3);

%%% Interpolating OISM data onto time grid
oism.datetime_interp = time_grid;
oism.temp1m_interp = interp1gap(datenum(oism.temp1m.datetime), oism.temp1m.temp, datenum(time_grid), datenum(max_gap));
oism.salt1m_interp = interp1gap(datenum(oism.salt1m.datetime), oism.salt1m.salt, datenum(time_grid), datenum(max_gap));
oism.temp7m_interp = interp1gap(datenum(oism.temp7m.datetime), oism.temp7m.temp, datenum(time_grid), datenum(max_gap));
oism.salt7m_interp = interp1gap(datenum(oism.salt7m.datetime), oism.salt7m.salt, datenum(time_grid), datenum(max_gap));
oism.pco2_7m_interp = interp1gap(datenum(oism.pco2_7m.datetime), oism.pco2_7m.pco2, datenum(time_grid), datenum(max_gap));
oism.do7m_interp = interp1gap(datenum(oism.do7m.datetime), oism.do7m.do, datenum(time_grid), datenum(max_gap));
oism.temp25m_interp = interp1gap(datenum(oism.temp25m.datetime), oism.temp25m.temp, datenum(time_grid), datenum(max_gap));
oism.salt25m_interp = interp1gap(datenum(oism.salt25m.datetime), oism.salt25m.salt, datenum(time_grid), datenum(max_gap));
oism.pco2_25m_interp = interp1gap(datenum(oism.pco2_25m.datetime), oism.pco2_25m.pco2, datenum(time_grid), datenum(max_gap));
oism.do25m_interp = interp1gap(datenum(oism.do25m.datetime), oism.do25m.do, datenum(time_grid), datenum(max_gap));

%%% Interpolating OOSM data onto time grid
ossm.datetime_interp = time_grid;
ossm.temp1m_interp = interp1gap(datenum(ossm.temp1m.datetime), ossm.temp1m.temp, datenum(time_grid), datenum(max_gap));
ossm.salt1m_interp = interp1gap(datenum(ossm.salt1m.datetime), ossm.salt1m.salt, datenum(time_grid), datenum(max_gap));
ossm.pco2_1m_interp = interp1gap(datenum(ossm.pco2_1m.datetime), ossm.pco2_1m.pco2, datenum(time_grid), datenum(max_gap));
ossm.temp7m_interp = interp1gap(datenum(ossm.temp7m.datetime), ossm.temp7m.temp, datenum(time_grid), datenum(max_gap));
ossm.salt7m_interp = interp1gap(datenum(ossm.salt7m.datetime), ossm.salt7m.salt, datenum(time_grid), datenum(max_gap));
ossm.do7m_interp = interp1gap(datenum(ossm.do7m.datetime), ossm.do7m.do, datenum(time_grid), datenum(max_gap));

%%% Interpolating OOSM data onto time grid
oosm.datetime_interp = time_grid;
oosm.temp1m_interp = interp1gap(datenum(oosm.temp1m.datetime), oosm.temp1m.temp, datenum(time_grid), datenum(max_gap));
oosm.salt1m_interp = interp1gap(datenum(oosm.salt1m.datetime), oosm.salt1m.salt, datenum(time_grid), datenum(max_gap));
oosm.pco2_1m_interp = interp1gap(datenum(oosm.pco2_1m.datetime), oosm.pco2_1m.pco2, datenum(time_grid), datenum(max_gap));
oosm.temp7m_interp = interp1gap(datenum(oosm.temp7m.datetime), oosm.temp7m.temp, datenum(time_grid), datenum(max_gap));
oosm.salt7m_interp = interp1gap(datenum(oosm.salt7m.datetime), oosm.salt7m.salt, datenum(time_grid), datenum(max_gap));
oosm.do7m_interp = interp1gap(datenum(oosm.do7m.datetime), oosm.do7m.do, datenum(time_grid), datenum(max_gap));

%%% Interpolating shelf wind data onto time grid
metero_shelf.datetime_interp = time_grid;
metero_shelf.wind_dir_interp = interp1gap(datenum(metero_shelf.datetime), metero_shelf.wind_dir, datenum(time_grid), datenum(max_gap));
metero_shelf.wind_spd_interp = interp1gap(datenum(metero_shelf.datetime), metero_shelf.wind_spd, datenum(time_grid), datenum(max_gap));

%%% Interpolating offshore wind data onto time grid
metero_offshore.datetime_interp = time_grid;
metero_offshore.wind_dir_interp = interp1gap(datenum(metero_offshore.datetime), metero_offshore.wind_dir, datenum(time_grid), datenum(max_gap));
metero_offshore.wind_spd_interp = interp1gap(datenum(metero_offshore.datetime), metero_offshore.wind_spd, datenum(time_grid), datenum(max_gap));

%%% Interpolating river discharge data
riverflow.datetime_interp = time_grid;
riverflow.flow_interp = interp1gap(datenum(riverflow.datetime), riverflow.flow, datenum(time_grid), datenum(max_gap));
    
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Creating a running mean for all data %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Setting the time step for the running mean
runmean_time_step = datenum(hours(1));

%%% Creating running mean for the OISM data on time grid
oism.datetime_runmean = time_grid;
oism.temp1m_runmean = runmean(time_grid, runmean_time_step, oism.temp1m.datetime, oism.temp1m.temp);
oism.salt1m_runmean = runmean(time_grid, runmean_time_step, oism.salt1m.datetime, oism.salt1m.salt);
oism.temp7m_runmean = runmean(time_grid, runmean_time_step, oism.temp7m.datetime, oism.temp7m.temp);
oism.salt7m_runmean = runmean(time_grid, runmean_time_step, oism.salt7m.datetime, oism.salt7m.salt);
oism.pco2_7m_runmean = runmean(time_grid, runmean_time_step, oism.pco2_7m.datetime, oism.pco2_7m.pco2);
oism.do7m_runmean = runmean(time_grid, runmean_time_step, oism.do7m.datetime, oism.do7m.do);
oism.temp25m_runmean = runmean(time_grid, runmean_time_step, oism.temp25m.datetime, oism.temp25m.temp);
oism.salt25m_runmean = runmean(time_grid, runmean_time_step, oism.salt25m.datetime, oism.salt25m.salt);
oism.pco2_25m_runmean = runmean(time_grid, runmean_time_step, oism.pco2_25m.datetime, oism.pco2_25m.pco2);
oism.do25m_runmean = runmean(time_grid, runmean_time_step, oism.do25m.datetime, oism.do25m.do);

%%% Creating running mean for the OSSM data on time grid
ossm.datetime_runmean = time_grid;
ossm.temp1m_runmean = runmean(time_grid, runmean_time_step, ossm.temp1m.datetime, ossm.temp1m.temp);
ossm.salt1m_runmean = runmean(time_grid, runmean_time_step, ossm.salt1m.datetime, ossm.salt1m.salt);
ossm.pco2_1m_runmean = runmean(time_grid, runmean_time_step, ossm.pco2_1m.datetime, ossm.pco2_1m.pco2);
ossm.temp7m_runmean = runmean(time_grid, runmean_time_step, ossm.temp7m.datetime, ossm.temp7m.temp);
ossm.salt7m_runmean = runmean(time_grid, runmean_time_step, ossm.salt7m.datetime, ossm.salt7m.salt);
ossm.do7m_runmean = runmean(time_grid, runmean_time_step, ossm.do7m.datetime, ossm.do7m.do);

%%% Creating running mean for the OOSM data on time grid
oosm.datetime_runmean = time_grid;
oosm.temp1m_runmean = runmean(time_grid, runmean_time_step, oosm.temp1m.datetime, oosm.temp1m.temp);
oosm.salt1m_runmean = runmean(time_grid, runmean_time_step, oosm.salt1m.datetime, oosm.salt1m.salt);
oosm.pco2_1m_runmean = runmean(time_grid, runmean_time_step, oosm.pco2_1m.datetime, oosm.pco2_1m.pco2);
oosm.temp7m_runmean = runmean(time_grid, runmean_time_step, oosm.temp7m.datetime, oosm.temp7m.temp);
oosm.salt7m_runmean = runmean(time_grid, runmean_time_step, oosm.salt7m.datetime, oosm.salt7m.salt);
oosm.do7m_runmean = runmean(time_grid, runmean_time_step, oosm.do7m.datetime, oosm.do7m.do);

%%% Creating running mean for shelf wind data on time grid
metero_shelf.datetime_runmean = time_grid;
metero_shelf.wind_dir_runmean = runmean(time_grid, runmean_time_step, metero_shelf.datetime, metero_shelf.wind_dir);
metero_shelf.wind_spd_runmean = runmean(time_grid, runmean_time_step, metero_shelf.datetime, metero_shelf.wind_spd);

%%% Creating running mean for offshore wind data on time grid
metero_offshore.datetime_runmean = time_grid;
metero_offshore.wind_dir_runmean = runmean(time_grid, runmean_time_step, metero_offshore.datetime, metero_offshore.wind_dir);
metero_offshore.wind_spd_runmean = runmean(time_grid, runmean_time_step, metero_offshore.datetime, metero_offshore.wind_spd);

%%% Creating running mean for river discharge data on time grid
riverflow.datetime_runmean = time_grid;
riverflow.flow_runmean = runmean(time_grid, runmean_time_step, riverflow.datetime, riverflow.flow);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Creating Figure %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [100 100 1000 1000])
sgtitle('Oregon Inshore Surface Mooring - 7m')

%%% Temperature Subplot
ax1 = subplot(311);
plot(oism.datetime_interp, oism.temp7m_interp)
ylabel('Temperature (degC)')

%%% Salinity Subplot
ax2 = subplot(312);
plot(oism.datetime_interp, oism.salt7m_interp);
ylabel('Salinity (psu)')

%%% pCO2 Subplot
ax3 = subplot(313);
plot(oism.datetime_interp, oism.pco2_7m_interp);
ylabel('pCO_2 (microatm)')

linkaxes([ax1 ax2 ax3], 'x')

clear ax1 ax2 ax3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [100 100 1000 1000])
sgtitle('Oregon Inshore Surface Mooring - 7m')

%%% Temperature Subplot
ax1 = subplot(311);
plot(oism.datetime_interp, oism.temp7m_interp)
ylabel('Temperature (degC)')
xlim([datetime(2014,04,17) datetime(2022,07,01)])

%%% Salinity Subplot
ax2 = subplot(312);
plot(oism.datetime_interp, oism.salt7m_interp);
ylabel('Salinity (psu)')
xlim([datetime(2014,04,17) datetime(2022,07,01)])

%%% DO Subplot
ax3 = subplot(313);
plot(oism.datetime_interp, oism.do7m_interp);
ylabel('Dissolved Oxygen (micromol/kg)')
xlim([datetime(2014,04,17) datetime(2022,07,01)])

linkaxes([ax1 ax2 ax3], 'x')

clear ax1 ax2 ax3
