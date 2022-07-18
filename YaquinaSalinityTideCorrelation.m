%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loading Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading Hatfield Time Series
load('YaquinaTS.mat');

%%% Loading Tide Data
load('HourlyTides.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Formatting Data on Time Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finding the first and last days of the Hatfield time series
firstday = min(yaquina_all.datenum);
lastday = max(yaquina_all.datenum);

%%% Creating a daily time grid for the Hatfield data
hourlytimegrid = floor(firstday):datenum(hours(1)):ceil(lastday);

%%% Creating a running mean for Hatfield salinity data
yaquina_all.datenum_1hrrunmean = hourlytimegrid;
yaquina_all.salt_1hrrunmean = runmean(hourlytimegrid, yaquina_all.datenum, yaquina_all.salt);

%%% Creating a running mean for tide data
tides_runmean.datenum_1hrrunmean = hourlytimegrid;
tides_runmean.MSL_1hrrunmean = runmean(hourlytimegrid, datenum(tides.datetime), tides.MSL);
