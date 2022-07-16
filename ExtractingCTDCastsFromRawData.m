%%% Importing CTD Cast data
castfull = readtable('065816_20201006_1510_data.txt');
name = '065816_20201006_1510_data';
i = 21;

%%% Plotting Pressure
figure('Renderer', 'painters', 'Position', [100 100 1000 1000])
plot(datenum(castfull.Time), castfull.SeaPressure)
datetick('x', 'hh:MM:ss');

%%% Pausing to zoom into cast
pause(20)

%%% Extracting start/end times of cast
[time, ~] = ginput(2);

%%% Extracting cast data within specified time interval
cast = castfull(datenum(castfull.Time) > time(1) & datenum(castfull.Time) < time(2), :);

%%% Assigning cast data to structure
casts_final(i).Name = name;
casts_final(i).Time = vertcat(cast.Time);
casts_final(i).Conductivity = vertcat(cast.Conductivity);
casts_final(i).Temperature = vertcat(cast.Temperature);
casts_final(i).Pressure = vertcat(cast.Pressure);
casts_final(i).Turbidity = vertcat(cast.Turbidity);
casts_final(i).SeaPressure = vertcat(cast.SeaPressure);
casts_final(i).Depth = vertcat(cast.Depth);
casts_final(i).Salinity = vertcat(cast.Salinity);
casts_final(i).SpecificConductivity = vertcat(cast.SpecificConductivity);
casts_final(i).SpeedOfSound = vertcat(cast.SpeedOfSound);
casts_final(i).DensityAnomaly = vertcat(cast.DensityAnomaly);

%%% Checking profile
plot(cast.Temperature, cast.SeaPressure)
set(gca, 'YDir', 'reverse')

clear castfull name time cast
