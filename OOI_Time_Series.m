%%% Loading data
temp_7m = table2struct(readtable('OOI_CE_OISM_NSIF_Temperature.csv'));
salt_7m = table2struct(readtable('OOI_CE_OISM_NSIF_Salinity.csv'));
temp_25m = table2struct(readtable('OOI_CE_OISM_SMFN_Temperature.csv'));
salt_25m = table2struct(readtable('OOI_CE_OISM_SMFN_Salinity.csv'));

%%% Filtering to data with the highest quality flag
temp_7m = temp_7m(vertcat(temp_7m.sea_water_temperature_qc_agg) == 1, :);
salt_7m = salt_7m(vertcat(salt_7m.sea_water_practical_salinity_qc_agg) == 1, :);
temp_25m = temp_25m(vertcat(temp_25m.sea_water_temperature_qc_agg) == 1, :);
salt_25m = salt_25m(vertcat(salt_25m.sea_water_practical_salinity_qc_agg) == 1, :);

%%% Converting time string to datetime
temp_7m_datetime = datetime(vertcat(temp_7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
salt_7m_datetime = datetime(vertcat(salt_7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
temp_25m_datetime = datetime(vertcat(temp_25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
salt_25m_datetime = datetime(vertcat(salt_25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

for i = 1:length(temp_7m_datetime)
    temp_7m(i).datetime = temp_7m_datetime(i);
end

for i = 1:length(salt_7m_datetime)
    salt_7m(i).datetime = salt_7m_datetime(i);
end
for i = 1:length(temp_25m_datetime)
    temp_25m(i).datetime = temp_25m_datetime(i);
end

for i = 1:length(salt_25m_datetime)
    salt_25m(i).datetime = salt_25m_datetime(i);
end

clear salt_7m_datetime temp_7m_datetime salt_25m_datetime temp_25m_datetime i 

%%
%%% Plotting temp/salt time series
figure()
sgtitle('Oregon Inshoor Surface Mooring - Near Surface Instrument Frame (7 m)');

subplot(211)
plot(vertcat(temp_7m.datetime), vertcat(temp_7m.sea_water_temperature))
ylabel('Temperature (degC)');

subplot(212)
plot(vertcat(salt_7m.datetime), vertcat(salt_7m.sea_water_practical_salinity));
ylabel('Salinity (g/kg)');

figure()
sgtitle('Oregon Inshoor Surface Mooring - Seafloor Multi-Function Node (25 m)');

subplot(211)
plot(vertcat(temp_25m.datetime), vertcat(temp_25m.sea_water_temperature))
ylabel('Temperature (degC)');

subplot(212)
plot(vertcat(salt_25m.datetime), vertcat(salt_25m.sea_water_practical_salinity));
ylabel('Salinity (g/kg)');
