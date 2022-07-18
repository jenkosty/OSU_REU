%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loading Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading Hatfield Time Series
load('YaquinaTS.mat');

%%% Loading OISM Time Series
load('OOI_CE_OISM.mat');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Finding Max Daily Salinity in Hatfield %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Finding the first and last days of the Hatfield time series
firstday = min(yaquina_all.datenum);
lastday = max(yaquina_all.datenum);

%%% Creating a daily time grid for the Hatfield data
dailytimegrid = floor(firstday):1:ceil(lastday);

clear firstday lastday

%%% Finding the max daily salinity
for i = length(dailytimegrid)-1:-1:1
    [maxdailysalinity(i).salt, idx] = max(yaquina_all.salt(yaquina_all.datenum > dailytimegrid(i) & yaquina_all.datenum < dailytimegrid(i+1)));
    dailydatetimes = yaquina_all.datenum(yaquina_all.datenum > dailytimegrid(i) & yaquina_all.datenum < dailytimegrid(i+1));
    maxdailysalinity(i).datenum = dailydatetimes(idx);
end

clear dailydatetimes idx i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Calculating High Tide Hatfield Temp/Salt %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(maxdailysalinity)
    
    %%% Extracting indices of measurements taken within 1 hour of highest
    %%% daily salinity
    idx = yaquina_all.datenum > maxdailysalinity(i).datenum-1 & yaquina_all.datenum < maxdailysalinity(i).datenum+1;
    
    %%% Calculating average salinity
    hatfieldHT.salt(i) = mean(yaquina_all.salt(idx), 'omitnan');
    
    %%% Calculating average temperature
    hatfieldHT.temp(i) = mean(yaquina_all.temp(idx), 'omitnan');
    
    %%% Extracting time of highest daily salinity
    hatfieldHT.datenum(i) = maxdailysalinity(i).datenum;
end

%%% Formatting as datetime
hatfieldHT.datetime = datetime(hatfieldHT.datenum, 'ConvertFrom', 'datenum');

clear idx i

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Interpolating Hatfield High Tide Data %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

noontimegrid = dailytimegrid+0.5;

hatfieldHT.datenum_interp = noontimegrid;
hatfieldHT.datetime_interp = datetime(hatfieldHT.datenum_interp, 'ConvertFrom', 'datenum');
hatfieldHT.salt_interp = interp1(hatfieldHT.datenum, hatfieldHT.salt, hatfieldHT.datenum_interp); 
hatfieldHT.temp_interp = interp1(hatfieldHT.datenum, hatfieldHT.temp, hatfieldHT.datenum_interp); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Calculating a Running Mean of OISM Data on Time Grid %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

oism.datetime_1dayrunmean = noontimegrid;

%%% 1m runmean for salinity and temperature
oism.salt1m_1dayrunmean = runmean(noontimegrid, datenum(oism.salt1m.datetime), oism.salt1m.salt);
oism.temp1m_1dayrunmean = runmean(noontimegrid, datenum(oism.temp1m.datetime), oism.temp1m.temp);

%%% 7m runmean for salinity and temperature
oism.salt7m_1dayrunmean = runmean(noontimegrid, datenum(oism.salt7m.datetime), oism.salt7m.salt);
oism.temp7m_1dayrunmean = runmean(noontimegrid, datenum(oism.temp7m.datetime), oism.temp7m.temp);

%%% Averaging 1m and 7m runmeans for salinity and temperature
oism.avgsalt_1dayrunmean = mean([oism.salt1m_1dayrunmean; oism.salt7m_1dayrunmean]);
oism.avgtemp_1dayrunmean = mean([oism.temp1m_1dayrunmean; oism.temp7m_1dayrunmean]);