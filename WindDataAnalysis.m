%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Loading Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Loading wind data
load('MeterologicalData.mat');

%%% Loading OISM data
load('OOI_CE_OISM.mat')

%%% Loading Yaquina River discharge data
load('YaquinaRiverDischarge.mat');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Decomposing Wind Data into U/V Components %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculating x and y components of the wind speed data
u = metero_shelf.wind_spd .* sin(metero_shelf.wind_dir .* pi/180);
v = metero_shelf.wind_spd .* cos(metero_shelf.wind_dir .* pi/180);

%%% Extracting summer values from Yaquina riverflow data %%%

%%% Identifying periods with low riverflow
summer_idx = riverflow.flow < 5;
summer_riverflow_flow_unfmt = riverflow.flow(summer_idx);
summer_riverflow_datetime_unfmt = riverflow.datetime(summer_idx);

summer_riverflow_flow = summer_riverflow_flow_unfmt;
summer_riverflow_datetime = summer_riverflow_datetime_unfmt;

%%% Adding between periods with large time gaps
for i = 1:length(summer_riverflow_datetime_unfmt)-1
    if diff([summer_riverflow_datetime_unfmt(i) summer_riverflow_datetime_unfmt(i+1)]) > days(1)
        summer_riverflow_datetime(end+1) = summer_riverflow_datetime_unfmt(i)+hours(12);
        summer_riverflow_flow(end+1) = NaN;
    end
end

%%% Sorting the NaNs to be in the appropriate vector location
[summer_riverflow_datetime, idy] = sort(summer_riverflow_datetime, 'ascend');
summer_riverflow_flow = summer_riverflow_flow(idy);

u_summer = u;
v_summer = v;

%%% Extracting summer u,v values for the wind data
for i = 1:length(metero_shelf.datetime)
    if min(abs(metero_shelf.datetime(i)-summer_riverflow_datetime)) > days(1)
        u_summer(i) = NaN;
        v_summer(i) = NaN;
    end
end

%%% Displaying u/v data
figure()
scatter(u_summer, v_summer)
title('Summer U/V Data');
        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Decomposing Wind Data into Principle/Secondary Components %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculating x and y components of the wind speed data
u = metero_shelf.wind_spd .* sin(metero_shelf.wind_dir .* pi/180);
v = metero_shelf.wind_spd .* cos(metero_shelf.wind_dir .* pi/180);

%%% Flipping the magnitude to adjust for meterological definiton of wind
%%% direction
u = -u;
v = -v;

%%% Finding the first and last days of the wind time series
first_day = min(metero_shelf.datenum);
last_day = max(metero_shelf.datenum);

%%% Creating a daily time grid for the wind data
dailytimegrid = floor(first_day):1:ceil(last_day);

%%% Creating a running mean for u/v data
u_rm = runmean(dailytimegrid, metero_shelf.datenum, u)';
v_rm = runmean(dailytimegrid, metero_shelf.datenum, v)';

%%% Identifying NaNs in data
nanIdx = isnan(u_rm) | isnan(v_rm); 

%%% Getting slope of line
a = v_rm(~nanIdx) \ u_rm(~nanIdx);

%%% Plotting data
figure()
hold on
scatter(v_rm(~nanIdx), u_rm(~nanIdx));
plot(v_rm, a.*v_rm, 'r');
grid on

%%% Extracting the principle wind component
b=90; b_rad=((b*pi)./180);
[THETA,R] = cart2pol(v_rm,a.*v_rm); %Convert to polar coordinates
THETA=THETA+b_rad; %Add a_rad to theta
[principle_dir_x,principle_dir_y] = pol2cart(THETA,R); %Convert back to Cartesian coordinates

%%% Extracting the secondary wind component
b=90; b_rad=((b*pi)./180);
[THETA,R] = cart2pol(principle_dir_x,principle_dir_y); %Convert to polar coordinates
THETA=THETA+b_rad; %Add a_rad to theta
[secondary_dir_x,secondary_dir_y] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
 
%%% Plotting data
figure()
hold on
scatter(u_rm(~nanIdx), v_rm(~nanIdx));
plot(principle_dir_x, principle_dir_y, 'r');
plot(secondary_dir_x, secondary_dir_y, 'g');
grid on

%%% Decomponsing the wind data into the principle and secondary components

principle_winds_prelim = NaN(size(v));
principle_winds_x = NaN(size(v));
principle_winds_y = NaN(size(v));
principle_winds = NaN(size(v));
secondary_winds_x = NaN(size(v));
secondary_winds_y = NaN(size(v));
secondary_winds = NaN(size(v));

for i = 1:length(v)
    
    %%% PRINCIPLE COMPONENT 
    %%% a dot b
    principle_winds_prelim(i) = u(i)*principle_dir_x(200) + v(i)*principle_dir_y(200);
    %%% Divide by b dot b
    principle_winds_prelim(i) = principle_winds_prelim(i) / (principle_dir_x(200)*principle_dir_x(200) + principle_dir_y(200)*principle_dir_y(200));
    %%% Multiply by b
    principle_winds_x(i) = principle_winds_prelim(i) * principle_dir_x(200);
    principle_winds_y(i) = principle_winds_prelim(i) * principle_dir_y(200);
    principle_winds(i) = sqrt(principle_winds_x(i)^2 + principle_winds_y(i)^2);
    if principle_winds_y(i) < 0
        principle_winds(i) = -1 * principle_winds(i);
    end
    
    %%% SECONDARY COMPONENT
    %%% a - proj(a)
    secondary_winds_x(i) = u(i) - principle_winds_x(i);
    secondary_winds_y(i) = v(i) - principle_winds_y(i);
    secondary_winds(i) = sqrt(secondary_winds_x(i)^2 + secondary_winds_y(i)^2);
    if secondary_winds_y(i) < 0
        secondary_winds(i) = -1 * secondary_winds(i);
    end
    
end

metero_shelf.u = u;
metero_shelf.v = v;
metero_shelf.principle_winds = principle_winds;
metero_shelf.secondary_winds = secondary_winds;

save('/Users/jenkosty/Research/OSU_REU/Processed_Data/MeterologicalData', 'metero_shelf', 'metero_offshore');
    
%% 

%%% Visualizing Decomposition Process

for i = 100
    figure()
    hold on
    
    xline(0, 'k');
    yline(0, 'k');
    scatter(u(i), v(i), 'k', 'MarkerFaceColor', 'k');
    plot(principle_dir_x, principle_dir_y, 'r');
    plot(secondary_dir_x, secondary_dir_y, 'g');
    scatter(principle_winds_x(i), principle_winds_y(i), 'r', 'MarkerFaceColor', 'r');
    scatter(secondary_winds_x(i), secondary_winds_y(i), 'g', 'MarkerFaceColor', 'g');
    grid on
    pause
    close all
end
    
    
%%

figure('Renderer', 'painters', 'Position', [100 100 1200 800])

ax1 = subplot(411);
hold on
area(metero_shelf.datetime, v, 'FaceColor', 'k', 'DisplayName', 'North-South Winds');
%plot(metero_shelf.datetime, principle_winds, 'k', 'DisplayName', 'Principle Winds')
yline(0, '--');
ylim([-20 20]);
ylabel('Wind Speed (m/s)');
title('Principle Winds');
legend()


ax2 = subplot(412);
hold on
area(metero_shelf.datetime, u, 'FaceColor', 'k','DisplayName', 'East-West Winds');
%plot(metero_shelf.datetime, secondary_winds, 'k', 'DisplayName', 'Secondary Winds');
yline(0, '--');
ylim([-20 20]);
ylabel('Wind Speed (m/s)');
title('Secondary Winds');
legend()

ax3 = subplot(413);
plot(oism.salt1m.datetime, oism.salt1m.salt, 'r')
ylabel('Salinity (psu)');
title('Inshore 1m Salinity');

ax4 = subplot(414);
hold on
plot(riverflow.datetime, riverflow.flow, 'b')
ylabel('River Discharge (m^3/s)');
title('Yaquina River Discharge');

linkaxes([ax1 ax2 ax3 ax4], 'x');
