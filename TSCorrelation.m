%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating a consistently spaced time grid
%time_grid = (datetime(2021,01,01):hours(2):datetime(2022,06,01))';
time_grid = yaquina_HT.datetime;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Loading and pre-processing Yaquina station data %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading station data
load('HFhightide.mat')

%%% Formatting in a structure
yaquina_HT.datetime = tm;
yaquina_HT.temp = T;
yaquina_HT.temp_std = Tsd;
yaquina_HT.salt = S;
yaquina_HT.salt_std = Ssd;

clear S Ssd T tm Tsd 

%%% Loading station data
load('HFsurface.mat')

%%% Formatting in a structure
yaquina_all.datetime = tm;
yaquina_all.temp = T;
yaquina_all.salt = S;
yaquina_all.cond = C;
yaquina_all.pres = P;

clear S T tm C P

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
oism.temp25m_runmean = runmean(time_grid, runmean_time_step, oism.temp25m.datetime, oism.temp25m.temp);
oism.salt25m_runmean = runmean(time_grid, runmean_time_step, oism.salt25m.datetime, oism.salt25m.salt);

%%% Creating running mean for the Yaquina Data
% yaquina_HT.datetime_runmean = time_grid;
% yaquina_HT.temp_runmean = runmean(time_grid, runmean_time_step, yaquina_HT.datetime, yaquina_HT.temp);
% yaquina_HT.salt_runmean = runmean(time_grid, runmean_time_step, yaquina_HT.datetime, yaquina_HT.salt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Conducting Correlation Analysis %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scatter(oism.temp1m_runmean,yaquina_HT.temp) 

hold on
plot(oism.datetime_runmean, oism.temp1m_runmean)
plot(yaquina_HT.datetime, yaquina_HT.temp)
hold off

[Rxy,mux,s2x,muy,s2y,k,Nk] = xcovar(oism.temp1m_runmean, yaquina_HT.temp, 10);

pearson_corr_coeff = Rxy./(s2x.*s2y);

[M,I] = max(pearson_corr_coeff);
