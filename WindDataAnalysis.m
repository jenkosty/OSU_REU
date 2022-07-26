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
%%%%%%%%%%%%%% Decomposing Wind Data into U/V Components %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculating x and y components of the wind speed data
u = metero_shelf.wind_spd .* sin(metero_shelf.wind_dir .* pi/180);
v = metero_shelf.wind_spd .* cos(metero_shelf.wind_dir .* pi/180);

%%% Finding the first and last days of the wind time series
first_day = min(metero_shelf.datenum);
last_day = max(metero_shelf.datenum);

%%% Creating a daily time grid for the wind data
dailytimegrid = floor(first_day):1:ceil(last_day);

%%% Creating a running mean for u/v data
u_rm = runmean(dailytimegrid, metero_shelf.datenum, u)';
v_rm = runmean(dailytimegrid, metero_shelf.datenum, v)';

nanIdx = isnan(u_rm) | isnan(v_rm); 
a = v_rm(~nanIdx) \ u_rm(~nanIdx);

p = polyfit(v_rm(~nanIdx), u_rm(~nanIdx), 1);
f = polyval(p, v_rm);

figure()
hold on
scatter(v_rm(~nanIdx), u_rm(~nanIdx));
plot(v_rm, a.*v_rm, 'r');
plot(v_rm, f, 'g');


b=90; b_rad=((b*pi)./180);
[THETA,R] = cart2pol(v_rm,a.*v_rm); %Convert to polar coordinates
THETA=THETA+b_rad; %Add a_rad to theta
[x,y] = pol2cart(THETA,R); %Convert back to Cartesian coordinates
 
 
figure()
hold on
scatter(u_rm(~nanIdx), v_rm(~nanIdx));
plot(x, y, 'r');

 
title('u/v data - daily running mean');

 



%%

figure('Renderer', 'painters', 'Position', [100 100 1200 800])

ax1 = subplot(311);
plot(metero_shelf.datetime, v, 'k')
ylim([-20 20]);
ylabel('Wind Speed (m/s)');
title('Wind Speed and Direction at the Shelf Station');

ax2 = subplot(312);
plot(oism.salt1m.datetime, oism.salt1m.salt, 'r')
ylabel('Salinity (psu)');
title('OISM 1m Salinity');

ax3 = subplot(313);
hold on
plot(riverflow.datetime, riverflow.flow, 'b')
plot(summer_riverflow_datetime, summer_riverflow_flow, 'r');
ylabel('River Discharge (m^3/s)');
title('Yaquina River Discharge');

linkaxes([ax1 ax2 ax3], 'x');

%%% Plot square axes of 



