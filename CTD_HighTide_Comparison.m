%%%%%%%%%%%%%%%%%
%%% Tide Data %%%
%%%%%%%%%%%%%%%%%

%%% Importing Yearly Tide Data
tides2020 = readtable('MTL_2020_HL.csv');
tides2021 = readtable('MTL_2021_HL.csv');
tides2022 = readtable('MTL_2022_HL.csv');

%%% Combining Tide Data into 1 table
tides = [tides2020; tides2021; tides2022];
clear tides2020 tides2021 tides2022

%%% Converting tide data to meters
tides.MTL = str2double(tides.Verified_ft_) / 3.281;

%%% Converting time data to datetime
times_string = string(tides.Date) + ' ' + string(tides.Time_LST_) + ':00';
tides.datetime = datetime(times_string, 'InputFormat', 'yyyy/MM/dd HH:mm:ss');
clear times_string

%%% Converting datetime data to UTC
tides.datetime = tides.datetime+hours(7);

%%%%%%%%%%%%%%%%
%%% CTD Data %%%
%%%%%%%%%%%%%%%%

%%% Importing CTD cast times
CTD_times = table2array(readtable('Hatfield_CTD_datetimes.txt'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Identifying CTD Casts Near High Tide %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Extracting high tide times
HT_times = tides.datetime(~isnan(tides.MTL) & tides.MTL>0);

%%% Finding the smallest duration between each CTD cast and high tide
for i = 1:length(CTD_times)
    time2HT(i,1) = min(abs(CTD_times(i) - HT_times));
end

%%% Saving CTD Data to a structure
CTDcasts = struct('datetime', CTD_times, 'time2hightide', time2HT);

clear HT_times time2HT HT_times CTD_times i

%%%%%%%%%%%%%%%%%%%%% 
%%% Plotting Data %%%
%%%%%%%%%%%%%%%%%%%%%
%%
figure('Renderer', 'painters', 'Position', [100 100 1500 500])
a = plot(tides.datetime(~isnan(tides.MTL)), tides.MTL(~isnan(tides.MTL)), 'k', 'DisplayName', 'H and L Tide'); % Plotting Tide Data
hold on

%%% Creating enteries for legend
b = xline(CTDcasts.datetime(1), 'm','DisplayName', 'CTD Casts');
c = xline(CTDcasts.datetime(1), 'c', 'DisplayName', 'CTD Casts Close to HT');

xline(CTDcasts.datetime, 'm') % Plotting CTD Cast Times
xline(CTDcasts.datetime(CTDcasts.time2hightide < hours(1)), 'c') % Plotting CTD Casts Close to High Tide

hold off
xlim([datetime(2020, 07, 01), datetime(2022, 07, 01)]);
ylabel('Mean Tide Level (m)')
legend([a,b,c])




