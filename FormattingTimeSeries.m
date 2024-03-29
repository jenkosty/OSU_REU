%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Creating a consistently spaced time grid
time_grid = (datetime(2014,04,01):days(1):datetime(2022,07,01))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Loading and pre-processing OOI data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OISM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Adding lat/lon data
oism.lat = 44.6598;
oism.lon = -124.095;

%%% Loading data
oism.temp1m = readtable('OOI_CE_OISM_SB_Temperature.csv');
oism.salt1m = readtable('OOI_CE_OISM_SB_Salinity.csv');
oism.pres1m = readtable('OOI_CE_OISM_SB_Pressure.csv');
oism.temp7m = readtable('OOI_CE_OISM_NSIF_Temperature.csv');
oism.salt7m = readtable('OOI_CE_OISM_NSIF_Salinity.csv');
oism.pres7m = readtable('OOI_CE_OISM_NSIF_Pressure.csv');
oism.density7m = readtable('OOI_CE_OISM_NSIF_Density.csv');
oism.pco2_7m = readtable('OOI_CE_OISM_NSIF_pCO2.csv');
oism.do7m = readtable('OOI_CE_OISM_NSIF_DO.csv');
oism.temp25m = readtable('OOI_CE_OISM_SMFN_Temperature.csv');
oism.salt25m = readtable('OOI_CE_OISM_SMFN_Salinity.csv');
oism.pres25m = readtable('OOI_CE_OISM_SMFN_Pressure.csv');
oism.density25m = readtable('OOI_CE_OISM_SMFN_Density.csv');
oism.pco2_25m = readtable('OOI_CE_OISM_SMFN_pCO2.csv');
oism.do25m = readtable('OOI_CE_OISM_SMFN_DO.csv');

%%% Filtering OISM to the highest quality control indicators
%%% Note: DO and pCO2 and density data do not have quality control indicators
oism.temp1m = oism.temp1m(oism.temp1m.sea_water_temperature_qc_agg == 1, :);
oism.salt1m = oism.salt1m(oism.salt1m.sea_water_practical_salinity_qc_agg == 1, :);
oism.pres1m = oism.pres1m(oism.pres1m.sea_water_pressure_qc_agg == 1, :);
oism.temp7m = oism.temp7m(oism.temp7m.sea_water_temperature_qc_agg == 1, :);
oism.salt7m = oism.salt7m(oism.salt7m.sea_water_practical_salinity_qc_agg == 1, :);
oism.pres7m = oism.pres7m(oism.pres7m.sea_water_pressure_qc_agg == 1, :);
oism.temp25m = oism.temp25m(oism.temp25m.sea_water_temperature_qc_agg == 1, :);
oism.salt25m = oism.salt25m(oism.salt25m.sea_water_practical_salinity_qc_agg == 1, :);
oism.pres25m = oism.pres25m(oism.pres25m.sea_water_pressure_qc_agg == 1, :);

%%% Converting time strings to datetime
oism.temp1m.datetime = datetime(vertcat(oism.temp1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt1m.datetime = datetime(vertcat(oism.salt1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.pres1m.datetime = datetime(vertcat(oism.pres1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.temp7m.datetime = datetime(vertcat(oism.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt7m.datetime = datetime(vertcat(oism.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.pres7m.datetime = datetime(vertcat(oism.pres7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.density7m.datetime = datetime(vertcat(oism.density7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.pco2_7m.datetime = datetime(vertcat(oism.pco2_7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.do7m.datetime = datetime(vertcat(oism.do7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.temp25m.datetime = datetime(vertcat(oism.temp25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt25m.datetime = datetime(vertcat(oism.salt25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.pres25m.datetime = datetime(vertcat(oism.pres25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.density25m.datetime = datetime(vertcat(oism.density25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.pco2_25m.datetime = datetime(vertcat(oism.pco2_25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.do25m.datetime = datetime(vertcat(oism.do25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
oism.temp1m = renamevars(oism.temp1m,'sea_water_temperature','temp');
oism.salt1m = renamevars(oism.salt1m,'sea_water_practical_salinity','salt');
oism.pres1m = renamevars(oism.pres1m, 'sea_water_pressure', 'pres');
oism.temp7m = renamevars(oism.temp7m,'sea_water_temperature','temp');
oism.salt7m = renamevars(oism.salt7m,'sea_water_practical_salinity','salt');
oism.pres7m = renamevars(oism.pres7m, 'sea_water_pressure', 'pres');
oism.density7m = renamevars(oism.density7m, 'sea_water_density', 'density');
oism.pco2_7m = renamevars(oism.pco2_7m,'partial_pressure_of_carbon_dioxide_in_sea_water','pco2');
oism.do7m = renamevars(oism.do7m,'moles_of_oxygen_per_unit_mass_in_sea_water','do');
oism.temp25m = renamevars(oism.temp25m,'sea_water_temperature','temp');
oism.salt25m = renamevars(oism.salt25m,'sea_water_practical_salinity','salt');
oism.pres25m = renamevars(oism.pres25m, 'sea_water_pressure', 'pres');
oism.density25m = renamevars(oism.density25m, 'sea_water_density', 'density');
oism.pco2_25m = renamevars(oism.pco2_25m,'partial_pressure_of_carbon_dioxide_in_sea_water','pco2');
oism.do25m = renamevars(oism.do25m,'moles_of_oxygen_per_unit_mass_in_sea_water','do');

%%% Mannually QCing salinity data
oism.salt25m.salt(oism.salt25m.datetime > datetime(2017, 11, 26) & oism.salt25m.datetime < datetime(2018, 1, 19)) = NaN; % in annotations
oism.salt25m.salt(oism.salt25m.datetime > datetime(2018, 11, 27) & oism.salt25m.datetime < datetime(2019, 4, 28)) = NaN;
oism.salt25m.salt(oism.salt25m.datetime > datetime(2019, 12, 13) & oism.salt25m.datetime < datetime(2020, 5, 30)) = NaN;
oism.salt25m.salt(oism.salt25m.datetime > datetime(2022, 1, 23) & oism.salt25m.datetime < datetime(2022, 4, 1)) = NaN;

%%% Mannually QCing temperature data
oism.temp25m.temp(oism.temp25m.datetime > datetime(2017, 11, 26) & oism.temp25m.datetime < datetime(2018, 1, 19)) = NaN; % in annotations

%%% Mannually QCing density data
oism.density25m.density(oism.density25m.datetime > datetime(2017, 11, 26) & oism.density25m.datetime < datetime(2018, 1, 19)) = NaN; % in annotations

%%% Mannually QCing pressure data
oism.pres25m.pres(oism.pres25m.datetime > datetime(2017, 11, 26) & oism.pres25m.datetime < datetime(2018, 1, 19)) = NaN; % in annotations

%%% Mannually QCing pCO2 data
oism.pco2_7m.pco2(oism.pco2_7m.pco2 > 2500) = NaN;
oism.pco2_25m.pco2(oism.pco2_25m.pco2 > 2500) = NaN;

%%% Removing pCO2 data flagged in OOI annotations
bad_data_7m = readtable('Bad_Data_OOI_CE_OISM_NSIF_pCO2.csv');
for i = 1:height(bad_data_7m)
    oism.pco2_7m.pco2(oism.pco2_7m.datetime > bad_data_7m.StartDate(i) & oism.pco2_7m.datetime < bad_data_7m.EndDate(i)) = NaN;
end

bad_data_25m = readtable('Bad_Data_OOI_CE_OISM_SMFN_pCO2.csv');
for i = 1:height(bad_data_25m)
    oism.pco2_25m.pco2(oism.pco2_25m.datetime > bad_data_25m.StartDate(i) & oism.pco2_25m.datetime < bad_data_25m.EndDate(i)) = NaN;
end

clear bad_data_7m bad_data_25m i

%%% Manually QCing DO data
oism.do7m.do(oism.do7m.do < 0) = NaN;
oism.do25m.do(oism.do25m.do < 0) = NaN;

%%% Removing DO data flagged in OOI annotations
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

%%% Removing Bad Data
oism.do7m_qc.do(ismember(oism.do7m_qc.datetime, oism.do7m_qc.fail_biofoul_datetime)) = NaN;

%%% Loading 25m QC DO Data
load('OOI_CE_OISM_SMFN_QC_DO.mat')

%%% Formatting
do25m_qc.do = dissolved_oxygen';
do25m_qc.datetime = datetime(time,'ConvertFrom','datenum')';
do25m_qc.fail_biofoul_datetime = datetime(fail_biofoul, 'ConvertFrom', 'datenum')';

%%% Saving to structure
oism.do25m_qc = do25m_qc;

%%% Removing Bad Data
oism.do7m_qc.do(ismember(oism.do25m_qc.datetime, oism.do25m_qc.fail_biofoul_datetime)) = NaN;

clear dissolved_oxygen fail_biofoul time do25m_qc

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/OOI_CE_OISM', 'oism')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OSSM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Adding lat/lon data
ossm.lat = 44.6393;
ossm.lon = -124.304;

%%% Loading data
ossm.temp1m = readtable('OOI_CE_OSSM_SB_Temperature.csv');
ossm.salt1m = readtable('OOI_CE_OSSM_SB_Salinity.csv');
ossm.pco2_1m = readtable('OOI_CE_OSSM_SB_pCO2.csv');
ossm.pco2air1m = readtable('OOI_CE_OSSM_SB_pCO2_Air.csv');
ossm.temp7m = readtable('OOI_CE_OSSM_NSIF_Temperature.csv');
ossm.salt7m = readtable('OOI_CE_OSSM_NSIF_Salinity.csv');
ossm.density7m = readtable('OOI_CE_OSSM_NSIF_Density.csv');
ossm.pres7m = readtable('OOI_CE_OSSM_NSIF_Pressure.csv');
ossm.do7m = readtable('OOI_CE_OSSM_NSIF_DO.csv');

%%% Filtering OSSM to the highest quality control indicators
%%% Note: DO and pCO2 data do not have quality control indicators
ossm.temp1m = ossm.temp1m(ossm.temp1m.sea_surface_temperature_qc_agg == 1, :);
ossm.salt1m = ossm.salt1m(ossm.salt1m.sea_water_practical_salinity_qc_agg == 1, :);
ossm.temp7m = ossm.temp7m(ossm.temp7m.sea_water_temperature_qc_agg == 1, :);
ossm.salt7m = ossm.salt7m(ossm.salt7m.sea_water_practical_salinity_qc_agg == 1, :);
ossm.pres7m = ossm.pres7m(ossm.pres7m.sea_water_pressure_qc_agg == 1, :);

%%% Converting time strings to datetime
ossm.temp1m.datetime = datetime(vertcat(ossm.temp1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.salt1m.datetime = datetime(vertcat(ossm.salt1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.pco2_1m.datetime = datetime(vertcat(ossm.pco2_1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.pco2air1m.datetime = datetime(vertcat(ossm.pco2air1m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.temp7m.datetime = datetime(vertcat(ossm.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.salt7m.datetime = datetime(vertcat(ossm.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.density7m.datetime = datetime(vertcat(ossm.density7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.pres7m.datetime = datetime(vertcat(ossm.pres7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
ossm.do7m.datetime = datetime(vertcat(ossm.do7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Renaming variables
ossm.temp1m = renamevars(ossm.temp1m,'sea_surface_temperature','temp');
ossm.salt1m = renamevars(ossm.salt1m,'sea_water_practical_salinity','salt');
ossm.pco2_1m = renamevars(ossm.pco2_1m,'surface_partial_pressure_of_carbon_dioxide_in_sea_water','pco2');
ossm.pco2air1m = renamevars(ossm.pco2air1m,'surface_partial_pressure_of_carbon_dioxide_in_air','pco2air');
ossm.temp7m = renamevars(ossm.temp7m,'sea_water_temperature','temp');
ossm.salt7m = renamevars(ossm.salt7m,'sea_water_practical_salinity','salt');
ossm.density7m = renamevars(ossm.density7m,'sea_water_density','density');
ossm.pres7m = renamevars(ossm.pres7m,'sea_water_pressure','pres');
ossm.do7m = renamevars(ossm.do7m,'moles_of_oxygen_per_unit_mass_in_sea_water','do');

%%% Manually QCing salinity data
ossm.salt1m.salt(ossm.salt1m.datetime > datetime(2022,7,23)) = NaN; % In annotations
ossm.salt1m.salt(ossm.salt1m.datetime > datetime(2021,12,16) & ossm.salt1m.datetime < datetime(2021,12,21)) = NaN; % In annotations
ossm.salt1m.salt(ossm.salt1m.datetime > datetime(2021,10,23) & ossm.salt1m.datetime < datetime(2021,10,26)) = NaN; % In annotations
ossm.salt1m.salt(ossm.salt1m.datetime > datetime(2021,9,19) & ossm.salt1m.datetime < datetime(2021,9,24)) = NaN; % In annotations
ossm.salt1m.salt(ossm.salt1m.datetime > datetime(2019,1,15) & ossm.salt1m.datetime < datetime(2019,4,20)) = NaN; % In annotations
ossm.salt1m.salt(ossm.salt1m.datetime > datetime(2017,4,20) & ossm.salt1m.datetime < datetime(2017,10,15)) = NaN; % In annotations

%%% Manually QCing temperature data
ossm.temp1m.temp(ossm.temp1m.datetime > datetime(2022,7,23)) = NaN; % In annotations
ossm.temp1m.temp(ossm.temp1m.datetime > datetime(2021,12,16) & ossm.temp1m.datetime < datetime(2021,12,21)) = NaN; % In annotations
ossm.temp1m.temp(ossm.temp1m.datetime > datetime(2021,10,23) & ossm.temp1m.datetime < datetime(2021,10,26)) = NaN; % In annotations
ossm.temp1m.temp(ossm.temp1m.datetime > datetime(2021,9,19) & ossm.temp1m.datetime < datetime(2021,9,24)) = NaN; % In annotations
ossm.temp1m.temp(ossm.temp1m.datetime > datetime(2019,1,15) & ossm.temp1m.datetime < datetime(2019,4,20)) = NaN; % In annotations
ossm.temp1m.temp(ossm.temp1m.datetime > datetime(2017,4,20) & ossm.temp1m.datetime < datetime(2017,10,15)) = NaN; % In annotations

%%% Manually QCing pCO2 data
ossm.pco2_1m.pco2(ossm.pco2_1m.pco2 > 2500) = NaN;

%%% Removing pCO2 data flagged in OOI annotations
bad_data_1m = readtable('Bad_Data_OOI_CE_OSSM_SB_pCO2.csv');
for i = 1:height(bad_data_1m)
    ossm.pco2_1m.pco2(ossm.pco2_1m.datetime > bad_data_1m.StartDate(i) & ossm.pco2_1m.datetime < bad_data_1m.EndDate(i)) = NaN;
end

%%% Removing atmospheric pCO2 data flagged in OOI annotations
bad_data_1m = readtable('Bad_Data_OOI_CE_OSSM_SB_pCO2_Air.csv');
for i = 1:height(bad_data_1m)
    ossm.pco2air1m.pco2(ossm.pco2air1m.datetime > bad_data_1m.StartDate(i) & ossm.pco2air1m.datetime < bad_data_1m.EndDate(i)) = NaN;
end

%%% Manually QCing DO data
ossm.do7m.do(ossm.do7m.do < 0) = NaN;

%%% Removing DO data flagged in OOI annotations
bad_data_7m = readtable('Bad_Data_OOI_CE_OSSM_NSIF_DO.csv');
for i = 1:height(bad_data_7m)
    ossm.do7m.do(ossm.do7m.datetime > bad_data_7m.StartDate(i) & ossm.do7m.datetime < bad_data_7m.EndDate(i)) = NaN;
end

clear bad_data_1m bad_data_7m i

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/OOI_CE_OSSM', 'ossm')
 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Loading OOSM data at different depths %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Adding lat/lon data
oosm.lat = 44.3811;
oosm.lon = -124.956;

%%% Loading data
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
oosm.pco2_1m.pco2(oosm.pco2_1m.pco2 > 2500) = NaN;

%%% Removing bad pCO2 data (according to OOI annotations)
bad_data_1m = readtable('Bad_Data_OOI_CE_OOSM_SB_pCO2.csv');
for i = 1:height(bad_data_1m)
    oosm.pco2_1m.pco2(oosm.pco2_1m.datetime > bad_data_1m.StartDate(i) & oosm.pco2_1m.datetime < bad_data_1m.EndDate(i)) = NaN;
end

%%% Filtering DO values to remove obviously bad data
oosm.do7m.do(oosm.do7m.do < 0) = NaN;

%%% Removing bad DO data (according to OOI annotations)
bad_data_7m = readtable('Bad_Data_OOI_CE_OOSM_NSIF_DO.csv');
for i = 1:height(bad_data_7m)
    oosm.do7m.do(oosm.do7m.datetime > bad_data_7m.StartDate(i) & oosm.do7m.datetime < bad_data_7m.EndDate(i)) = NaN;
end

clear bad_data_1m bad_data_7m i

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/OOI_CE_OOSM', 'oosm')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Loading and pre-processing meterological data %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%
%%% Shelf Data %%%
%%%%%%%%%%%%%%%%%%

%%% Adding lat/lon data
metero_shelf.lat = 44.639;
metero_shelf.lon = -124.304;

%%% Loading meterological data from different years
metero_unfmt.data2016 = readtable('Shelf Wind Data/2016meterologicaldata.txt');
metero_unfmt.data2017 = readtable('Shelf Wind Data/2017meterologicaldata.txt');
metero_unfmt.data2018 = readtable('Shelf Wind Data/2018meterologicaldata.txt');
metero_unfmt.data2019 = readtable('Shelf Wind Data/2019meterologicaldata.txt');
metero_unfmt.data2020 = readtable('Shelf Wind Data/2020meterologicaldata.txt');
metero_unfmt.data2021 = readtable('Shelf Wind Data/2021meterologicaldata.txt');
metero_unfmt.dataJan2022 = readtable('Shelf Wind Data/Jan2022meterologicaldata.txt');
metero_unfmt.dataFeb2022 = readtable('Shelf Wind Data/Feb2022meterologicaldata.txt');
metero_unfmt.dataMar2022 = readtable('Shelf Wind Data/Mar2022meterologicaldata.txt');
metero_unfmt.dataApr2022 = readtable('Shelf Wind Data/Apr2022meterologicaldata.txt');
metero_unfmt.dataMay2022 = readtable('Shelf Wind Data/May2022meterologicaldata.txt');
metero_unfmt.dataJune2022 = readtable('Shelf Wind Data/June2022meterologicaldata.txt');
metero_unfmt.dataJulyAug2022 = readtable('Shelf Wind Data/JulyAug2022meterologicaldata.txt');

%%% Removing June data from July/Augest section
metero_unfmt.dataJulyAug2022 = metero_unfmt.dataJulyAug2022(metero_unfmt.dataJulyAug2022.MM > 6, :);

%%% Extracting wind direction data
metero_shelf.wind_dir = table2array(vertcat(metero_unfmt.data2016(:,6),metero_unfmt.data2017(:,6),metero_unfmt.data2018(:,6),...
    metero_unfmt.data2019(:,6), metero_unfmt.data2020(:,6), metero_unfmt.data2021(:,6), metero_unfmt.dataJan2022(:,6),...
    metero_unfmt.dataFeb2022(:,6),metero_unfmt.dataMar2022(:,6), metero_unfmt.dataApr2022(:,6), metero_unfmt.dataMay2022(:,6),...
    metero_unfmt.dataJune2022(:,6), metero_unfmt.dataJulyAug2022(:,6)));

%%% Extracting wind speed data
metero_shelf.wind_spd = table2array(vertcat(metero_unfmt.data2016(:,7),metero_unfmt.data2017(:,7),metero_unfmt.data2018(:,7),...
    metero_unfmt.data2019(:,7), metero_unfmt.data2020(:,7), metero_unfmt.data2021(:,7), metero_unfmt.dataJan2022(:,7),...
    metero_unfmt.dataFeb2022(:,7), metero_unfmt.dataMar2022(:,7), metero_unfmt.dataApr2022(:,7), metero_unfmt.dataMay2022(:,7),...
    metero_unfmt.dataJune2022(:,7), metero_unfmt.dataJulyAug2022(:,7)));

%%% Extracting time of measurement data
yrs = table2array(vertcat(metero_unfmt.data2016(:,1),metero_unfmt.data2017(:,1),metero_unfmt.data2018(:,1),...
    metero_unfmt.data2019(:,1), metero_unfmt.data2020(:,1), metero_unfmt.data2021(:,1), metero_unfmt.dataJan2022(:,1),...
    metero_unfmt.dataFeb2022(:,1),metero_unfmt.dataMar2022(:,1), metero_unfmt.dataApr2022(:,1), metero_unfmt.dataMay2022(:,1),...
    metero_unfmt.dataJune2022(:,1), metero_unfmt.dataJulyAug2022(:,1)));

mths = table2array(vertcat(metero_unfmt.data2016(:,2),metero_unfmt.data2017(:,2),metero_unfmt.data2018(:,2),...
    metero_unfmt.data2019(:,2), metero_unfmt.data2020(:,2), metero_unfmt.data2021(:,2), metero_unfmt.dataJan2022(:,2),...
    metero_unfmt.dataFeb2022(:,2),metero_unfmt.dataMar2022(:,2), metero_unfmt.dataApr2022(:,2), metero_unfmt.dataMay2022(:,2),...
    metero_unfmt.dataJune2022(:,2), metero_unfmt.dataJulyAug2022(:,2)));

days = table2array(vertcat(metero_unfmt.data2016(:,3),metero_unfmt.data2017(:,3),metero_unfmt.data2018(:,3),...
    metero_unfmt.data2019(:,3), metero_unfmt.data2020(:,3), metero_unfmt.data2021(:,3), metero_unfmt.dataJan2022(:,3),...
    metero_unfmt.dataFeb2022(:,3),metero_unfmt.dataMar2022(:,3), metero_unfmt.dataApr2022(:,3), metero_unfmt.dataMay2022(:,3),...
    metero_unfmt.dataJune2022(:,3), metero_unfmt.dataJulyAug2022(:,3)));

hrs = table2array(vertcat(metero_unfmt.data2016(:,4),metero_unfmt.data2017(:,4),metero_unfmt.data2018(:,4),...
    metero_unfmt.data2019(:,4), metero_unfmt.data2020(:,4), metero_unfmt.data2021(:,4), metero_unfmt.dataJan2022(:,4),...
    metero_unfmt.dataFeb2022(:,4),metero_unfmt.dataMar2022(:,4), metero_unfmt.dataApr2022(:,4), metero_unfmt.dataMay2022(:,4),...
    metero_unfmt.dataJune2022(:,4), metero_unfmt.dataJulyAug2022(:,4)));

mnts = table2array(vertcat(metero_unfmt.data2016(:,5),metero_unfmt.data2017(:,5),metero_unfmt.data2018(:,5),...
    metero_unfmt.data2019(:,5), metero_unfmt.data2020(:,5), metero_unfmt.data2021(:,5), metero_unfmt.dataJan2022(:,5),...
    metero_unfmt.dataFeb2022(:,5),metero_unfmt.dataMar2022(:,5), metero_unfmt.dataApr2022(:,5), metero_unfmt.dataMay2022(:,5),...
    metero_unfmt.dataJune2022(:,5), metero_unfmt.dataJulyAug2022(:,5)));

snds = zeros(size(yrs)); %%% Setting seconds to 0

%%% Converting time of measurement data to datetime
metero_shelf.datetime = datetime(yrs,mths,days,hrs,mnts,snds);
metero_shelf.datenum = datenum(metero_shelf.datetime);

%%% Removing missing data
bad_data = find(metero_shelf.wind_spd > 90 | metero_shelf.wind_dir > 900);
metero_shelf.wind_dir(bad_data) = [];
metero_shelf.wind_spd(bad_data) = [];
metero_shelf.datetime(bad_data) = [];
metero_shelf.datenum(bad_data) = [];

clear yrs mths days hrs mnts snds metero_unfmt bad_data

%%%%%%%%%%%%%%%%%%%%%
%%% Offshore Data %%%
%%%%%%%%%%%%%%%%%%%%%

%%% Adding lat/lon data
metero_offshore.lat = 44.669;
metero_offshore.lon = -124.546;

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
metero_offshore.datetime = datetime(yrs,mths,days,hrs,mnts,snds);
metero_offshore.datenum = datenum(metero_offshore.datetime);

%%% Removing missing data
bad_data = find(metero_offshore.wind_spd > 90 | metero_offshore.wind_dir > 900);
metero_offshore.wind_dir(bad_data) = [];
metero_offshore.wind_spd(bad_data) = [];
metero_offshore.datetime(bad_data) = [];
metero_offshore.datenum(bad_data) = [];

clear yrs mths days hrs mnts snds metero_unfmt bad_data

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/MeterologicalData', 'metero_shelf', 'metero_offshore');
%%
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
riverflow.datenum = datenum(riverflow.datetime);

%%% Adding lat/lon data
riverflow.lat = 44.657397;
riverflow.lon = -123.83875;

clear riverflow_unfmt

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/YaquinaRiverDischarge', 'riverflow');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Loading and pre-processing Yaquina station data %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading station data
load('HFhightide.mat')

%%% Adding lat/lon data
yaquina_HT.lat =  44.624504;
yaquina_HT.lon = -124.043125;

%%% Formatting in a structure
yaquina_HT.datetime = tm;
yaquina_HT.datenum = datenum(yaquina_HT.datetime);
yaquina_HT.temp = T;
yaquina_HT.temp_std = Tsd;
yaquina_HT.salt = S;
yaquina_HT.salt_std = Ssd;

clear S Ssd T tm Tsd 

%%% Loading station data
load('HFsurface.mat')

%%% Adding lat/lon data
yaquina_all.lat =  44.624504;
yaquina_all.lon = -124.043125;

%%% Formatting in a structure
yaquina_all.datetime = tm;
yaquina_all.datenum = datenum(yaquina_all.datetime);
yaquina_all.temp = T;
yaquina_all.salt = S;
yaquina_all.cond = C;
yaquina_all.pres = P;

clear S T tm C P

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/YaquinaTS', 'yaquina_HT', 'yaquina_all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Loading and pre-processing Dock 5 data %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading station data
load('DOCK5_proc01.mat')

dock5.datenum = DND5;
dock5.datetime = datetime(dock5.datenum, 'ConvertFrom', 'datenum');
dock5.temp3ft = TD5(:,1);
dock5.salt3ft = SD5(:,1);
dock5.temp7ft = TD5(:,2);
dock5.salt7ft = SD5(:,2);
dock5.temp11ft = TD5(:,3);
dock5.salt11ft = SD5(:,3);

clear DND5 SD5 TD5

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/Dock5', 'dock5');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Loading and pre-processing tide data %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%% Importing Yearly Tide Data
tides2020 = readtable('MSL_2020_HL.csv');
tides2021 = readtable('MSL_2021_HL.csv');
tides2022 = readtable('MSL_2022_HL.csv');

%%% Combining Tide Data into 1 table
tides_unfmt = [tides2020; tides2021; tides2022];
clear tides2020 tides2021 tides2022

%%% Extracting verified tide data
tides_unfmt.MSL = tides_unfmt.Verified_m_;

%%% Using predicted data if verified is not available
%tides_unfmt.MSL(isnan(tides_unfmt.MSL)) = tides_unfmt.Predicted_m_(isnan(tides_unfmt.MSL));

%%% Converting time data to datetime
times_string = string(tides_unfmt.Date) + ' ' + string(tides_unfmt.Time_GMT_) + ':00';
tides_unfmt.datetime = datetime(times_string, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
tides_unfmt.datenum = datenum(tides_unfmt.datetime);
clear times_string

%%% Converting to structure
tides.datetime = tides_unfmt.datetime;
tides.datenum = tides_unfmt.datenum;
tides.MSL = tides_unfmt.MSL;

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/HLTides', 'tides')

clear tides_unfmt

%%
%%% Importing Yearly Tide Data
tides2020 = readtable('MSL_2020_Hourly.csv');
tides2021 = readtable('MSL_2021_Hourly.csv');
tides2022 = readtable('MSL_2022_Hourly.csv');

%%% Combining Tide Data into 1 table
tides_unfmt = [tides2020; tides2021; tides2022];
clear tides2020 tides2021 tides2022

%%% Extracting verified tide data
tides_unfmt.MSL = tides_unfmt.Verified_m_;

%%% Using predicted data if verified is not available
tides_unfmt.MSL(isnan(tides_unfmt.MSL)) = tides_unfmt.Predicted_m_(isnan(tides_unfmt.MSL));

%%% Converting time data to datetime
times_string = string(tides_unfmt.Date) + ' ' + string(tides_unfmt.Time_GMT_) + ':00';
tides_unfmt.datetime = datetime(times_string, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
tides_unfmt.datenum = datenum(tides_unfmt.datetime);
clear times_string

%%% Converting to structure
tides.datetime = tides_unfmt.datetime;
tides.datenum = tides_unfmt.datenum;
tides.MSL = tides_unfmt.MSL;

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/HourlyTides', 'tides')

clear tides_unfmt

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
runmean_time_step = datenum(hours(12));

%%% Creating running mean for the OISM data on time grid
oism.datetime_1dayrunmean = time_grid;
oism.temp1m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.temp1m.datetime), oism.temp1m.temp);
oism.salt1m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.salt1m.datetime), oism.salt1m.salt);
oism.temp7m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.temp7m.datetime), oism.temp7m.temp);
oism.salt7m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.salt7m.datetime), oism.salt7m.salt);
oism.pco2_7m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.pco2_7m.datetime), oism.pco2_7m.pco2);
oism.do7m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.do7m.datetime), oism.do7m.do);
oism.temp25m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.temp25m.datetime), oism.temp25m.temp);
oism.salt25m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.salt25m.datetime), oism.salt25m.salt);
oism.pco2_25m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.pco2_25m.datetime), oism.pco2_25m.pco2);
oism.do25m_1dayrunmean = runmean(datenum(time_grid), datenum(oism.do25m.datetime), oism.do25m.do);

%%% Creating running mean for the OSSM data on time grid
ossm.datetime_1dayrunmean = time_grid;
ossm.temp1m_1dayrunmean = runmean(datenum(time_grid), datenum(ossm.temp1m.datetime), ossm.temp1m.temp);
ossm.salt1m_1dayrunmean = runmean(datenum(time_grid), datenum(ossm.salt1m.datetime), ossm.salt1m.salt);
ossm.pco2_1m_1dayrunmean = runmean(datenum(time_grid), datenum(ossm.pco2_1m.datetime), ossm.pco2_1m.pco2);
ossm.temp7m_1dayrunmean = runmean(datenum(time_grid), datenum(ossm.temp7m.datetime), ossm.temp7m.temp);
ossm.salt7m_1dayrunmean = runmean(datenum(time_grid), datenum(ossm.salt7m.datetime), ossm.salt7m.salt);
ossm.do7m_1dayrunmean = runmean(datenum(time_grid), datenum(ossm.do7m.datetime), ossm.do7m.do);

