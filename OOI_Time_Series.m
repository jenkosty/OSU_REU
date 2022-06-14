%%% Loading data
oism.temp7m = readtable('OOI_CE_OISM_NSIF_Temperature.csv');
oism.salt7m = readtable('OOI_CE_OISM_NSIF_Salinity.csv');
oism.temp25m = readtable('OOI_CE_OISM_SMFN_Temperature.csv');
oism.salt25m = readtable('OOI_CE_OISM_SMFN_Salinity.csv');

%%% Filtering Data
oism.temp7m = oism.temp7m(oism.temp7m.sea_water_temperature_qc_agg == 1, :);
oism.salt7m = oism.salt7m(oism.salt7m.sea_water_practical_salinity_qc_agg == 1, :);
oism.temp25m = oism.temp25m(oism.temp25m.sea_water_temperature_qc_agg == 1, :);
oism.salt25m = oism.salt25m(oism.salt25m.sea_water_practical_salinity_qc_agg == 1, :);

%%% Creating a consistently spaced time grid
time_grid = (datetime(2014,04,01):minutes(30):datetime(2022,06,01))';

%%% Converting time string to datetime
oism.temp7m.datetime = datetime(vertcat(oism.temp7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt7m.datetime = datetime(vertcat(oism.salt7m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.temp25m.datetime = datetime(vertcat(oism.temp25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');
oism.salt25m.datetime = datetime(vertcat(oism.salt25m.time), 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss''Z');

%%% Interpolating onto time grid
oism_interp.time = time_grid;
oism_interp.temp7m = interp1(datenum(oism.temp7m.datetime), oism.temp7m.sea_water_temperature, datenum(time_grid));
oism_interp.salt7m = interp1(datenum(oism.salt7m.datetime), oism.salt7m.sea_water_practical_salinity, datenum(time_grid));
oism_interp.temp25m = interp1(datenum(oism.temp25m.datetime), oism.temp25m.sea_water_temperature, datenum(time_grid));
oism_interp.salt25m = interp1(datenum(oism.salt25m.datetime), oism.salt25m.sea_water_practical_salinity, datenum(time_grid));

%%% Removing 
for i = 1:length(oism_interp.time)
    if min(abs(datenum(oism_interp.time(i))-datenum(oism.temp7m.datetime))) > datenum(hours(1))
        oism_interp.temp7m(i) = NaN; 
    end
    if min(abs(datenum(oism_interp.time(i))-datenum(oism.salt7m.datetime))) > datenum(hours(1))
        oism_interp.salt7m(i) = NaN; 
    end
    if min(abs(datenum(oism_interp.time(i))-datenum(oism.temp25m.datetime))) > datenum(hours(1))
        oism_interp.temp25m(i) = NaN; 
    end
    if min(abs(datenum(oism_interp.time(i))-datenum(oism.salt25m.datetime))) > datenum(hours(1))
        oism_interp.salt25m(i) = NaN; 
    end
end

%%
%%% Plotting temp/salt time series
figure()
sgtitle('Oregon Inshoor Surface Mooring - Near Surface Instrument Frame (7 m)');

subplot(211)
plot(oism_interp.time, oism_interp.temp7m)
ylabel('Temperature (degC)');

subplot(212)
plot(oism_interp.time, oism_interp.salt7m)
ylabel('Salinity (g/kg)');

figure()
sgtitle('Oregon Inshoor Surface Mooring - Seafloor Multi-Function Node (25 m)');

subplot(211)
plot(oism_interp.time, oism_interp.temp25m);
ylabel('Temperature (degC)');

subplot(212)
plot(oism_interp.time, oism_interp.salt25m);
ylabel('Salinity (g/kg)');


