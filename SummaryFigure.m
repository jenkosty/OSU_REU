=%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Loading and pre-processing data %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading OISM data
load('OOI_CE_OISM.mat');

%%% Loading Yaquina time series
load('YaquinaTS.mat');

%%% Loading CTD Cast data
load('HatfieldCTDcasts.mat')

%%% Extracting high tide casts
casts_HT = casts_final(vertcat(casts_final.TimeToClosestHighTide) < hours(1));

%%% Loading tide data
load('HourlyTides.mat');

%%% Loading Yaquina River discharge data
load('YaquinaRiverDischarge.mat');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Creating Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Renderer', 'painters', 'Position', [100 100 1200 1000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Temperature subplot %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax1 = subplot(411);
hold on
a = plot(yaquina_all.datetime, yaquina_all.temp, 'Color', [.7 .7 .7], 'DisplayName', 'Hatfield Full Resolution');
b = plot(oism.temp1m.datetime, oism.temp1m.temp, 'Color', [68 68 253]./255, 'DisplayName', 'OISM (1m)');
c = plot(yaquina_HT.datetime, yaquina_HT.temp, 'k', 'DisplayName', 'Hatfield High Tide');
for i = 1:length(casts_HT)   

    %%% Plotting high tide cast
    d = plot([mean(casts_HT(i).Time) mean(casts_HT(i).Time)], ylim, '--k', 'DisplayName', 'High Tide CTD Casts');
    uistack(d, "bottom");
    
    %%% Getting indices for top/bottom of cast
    cast_top = casts_HT(i).Depth > min(casts_HT(i).Depth) & casts_HT(i).Depth < 1;
    cast_bottom = casts_HT(i).Depth < max(casts_HT(i).Depth) & casts_HT(i).Depth > max(casts_HT(i).Depth)-1;
    
    %%% Calculating mean temperature for top/bottom of cast
    cast_top_temp = mean(casts_HT(i).Temperature(cast_top));
    cast_bottom_temp = mean(casts_HT(i).Temperature(cast_bottom));
    
    %%% Plotting mean top/bottom temperatures of cast
    e = plot(mean(casts_HT(i).Time), cast_top_temp, 'ok', 'MarkerFaceColor', 'm','MarkerSize', 9, 'DisplayName', 'Top Temperature');
    f = plot(mean(casts_HT(i).Time), cast_bottom_temp, 'ok', 'MarkerFaceColor', 'g','MarkerSize', 9, 'DisplayName', 'Bottom Temperature');
    
end
hold off
xlim([datetime(2020,12,1), datetime(2022,6,1)])
ylabel('Temperature (degC)')
legend([a b c d e f], 'Location', 'southwest', 'Orientation', 'horizontal');

clear cast_top_temp cast_bottom_temp cast_top cast_bottom a b c d e f

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Salinity subplot %%%
%%%%%%%%%%%%%%%%%%%%%%%%

ax2 = subplot(412);
hold on
a = plot(yaquina_all.datetime, yaquina_all.salt, 'Color', [.7 .7 .7], 'DisplayName', 'Hatfield Full Resolution');
b = plot(oism.salt1m.datetime, oism.salt1m.salt, 'Color', [68 68 253]./255, 'DisplayName', 'OISM (1m)');
c = plot(yaquina_HT.datetime, yaquina_HT.salt, 'k', 'DisplayName', 'Hatfield High Tide');
for i = 1:length(casts_HT)  
    
    %%% Plotting high tide cast
    d = plot([mean(casts_HT(i).Time) mean(casts_HT(i).Time)], ylim, '--k', 'DisplayName', 'High Tide CTD Casts');
    uistack(d, "bottom");
    
    %%% Getting indices for top/bottom of cast
    cast_top = casts_HT(i).Depth > min(casts_HT(i).Depth) & casts_HT(i).Depth < 1;
    cast_bottom = casts_HT(i).Depth < max(casts_HT(i).Depth) & casts_HT(i).Depth > max(casts_HT(i).Depth)-1;
    
    %%% Calculating mean salinities for top/bottom of cast
    cast_top_salt = mean(casts_HT(i).Salinity(cast_top));
    cast_bottom_salt = mean(casts_HT(i).Salinity(cast_bottom));
    
    %%% Plotting mean top/bottom salinities of cast
    e = plot(mean(casts_HT(i).Time), cast_top_salt, 'ok', 'MarkerFaceColor', 'm','MarkerSize', 9, 'DisplayName', 'Top Salinity');
    f = plot(mean(casts_HT(i).Time), cast_bottom_salt, 'ok', 'MarkerFaceColor', 'g','MarkerSize', 9, 'DisplayName', 'Bottom Salinity');
    
end
hold off
xlim([datetime(2020,12,1), datetime(2022,6,1)])
ylabel('Salinity (psu)');
legend([a b c d e f], 'Location', 'southwest', 'Orientation', 'horizontal');

clear cast_top_salt cast_bottom_salt cast_top cast_bottom i a b c d e f

%%%%%%%%%%%%%%%%%%%%
%%% Tide Subplot %%%
%%%%%%%%%%%%%%%%%%%%

ax3 = subplot(413);
plot(tides.datetime, tides.MSL, 'k');
xlim([datetime(2020,12,1), datetime(2022,6,1)])
ylabel('Mean Sea Level (m)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Yaquina River discharge %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax4 = subplot(414);
plot(riverflow.datetime, riverflow.flow, 'k');
xlim([datetime(2020,12,1), datetime(2022,6,1)])
ylabel('Yaquina River Discharge (m^3/s)');

%%% Linking axes
linkaxes([ax1 ax2 ax3 ax4], 'x');
clear ax1 ax2 ax3 ax4
