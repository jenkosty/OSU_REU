%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Loading and pre-processing OISM data %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading OISM data
load('OOI_CE_OISM.mat');

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
%%%%%%%%%%%%%%%% Loading and Pre-Processing CTD Cast Data %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading cast data
load('HatfieldCTDcasts.mat')

%%% Extracting high tide casts
casts_HT = casts_final(vertcat(casts_final.TimeToClosestHighTide) < hours(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Creating Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
figure('Renderer', 'painters', 'Position', [100 100 1200 700])
sgtitle('CTD Casts Taken Within 2 Hours of High Tide')

%%% Temperature subplot
subplot(211)
hold on
a = plot(yaquina_all.datetime, yaquina_all.temp, 'Color', [.7 .7 .7], 'DisplayName', 'Hatfield - Full Resolution');
b = plot(oism.temp1m.datetime, oism.temp1m.temp, 'Color', [68 68 253]./255, 'DisplayName', 'OISM - 1m');
c = plot(yaquina_HT.datetime, yaquina_HT.temp, 'k', 'DisplayName', 'Hatfield - High Tide');
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
legend([a b c d e f]);

clear cast_top_temp cast_bottom_temp cast_top cast_bottom a b

%%% Salinity subplot
subplot(212)
hold on
a = plot(yaquina_all.datetime, yaquina_all.salt, 'Color', [.7 .7 .7], 'DisplayName', 'Hatfield - Full Resolution');
b = plot(oism.salt1m.datetime, oism.salt1m.salt, 'Color', [68 68 253]./255, 'DisplayName', 'OISM - 1m');
c = plot(yaquina_HT.datetime, yaquina_HT.salt, 'k', 'DisplayName', 'Hatfield - High Tide');
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
legend([a b c d e f], 'Location', 'southeast');

clear cast_top_salt cast_bottom_salt cast_top cast_bottom i a b c d e
