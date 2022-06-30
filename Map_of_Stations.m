%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Creating Map to Display Monitoring Stations %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Yaquina Station Coordinates
yaquina_station.lat =  44.624504;
yaquina_station.lon = -124.043125;

%%% OISM Coordinates
oism.lat = 44.6598;
oism.lon = -124.095;

%%% OSSM Coordinates
ossm.lat = 44.6393;
ossm.lon =  -124.304;

%%% OOSM Coordinates
oosm.lat = 44.3811;
oosm.lon = -124.956;

%%% Meterological Station Coordiantes
metero_shelf.lat = 44.639;
metero_shelf.lon = -124.304;

metero_offshore.lat = 44.669;
metero_offshore.lon = -124.546;

%%% River Discharge Station Coordinates
riverflow.lat = 44.657397;
riverflow.lon = -123.83875;

%%%

figure1 = figure('Renderer', 'painters', 'Position', [100 100 900 450]);
c = geoscatter([metero_shelf.lat metero_offshore.lat], [metero_shelf.lon metero_offshore.lon], 300, 'sk', 'filled', 'DisplayName', 'NDBC Wind Gauges');
hold on
d = geoscatter(riverflow.lat, riverflow.lon, 100, '^r', 'filled', 'DisplayName', 'Riverflow Gauge');
a = geoscatter([oism.lat ossm.lat oosm.lat],[oism.lon ossm.lon oosm.lon], 100, 'g', 'filled', 'DisplayName', 'OOI Surface Moorings');
text(oism.lat+0.005, oism.lon+0.015, 'Inshore', 'FontSize', 15);
text(ossm.lat+0.005, ossm.lon+0.015, 'Shelf', 'FontSize', 15)
text(oosm.lat+0.005, oosm.lon+0.015, 'Offshore', 'FontSize', 15)
b = geoscatter(yaquina_station.lat,yaquina_station.lon, 100, 'b', 'filled','DisplayName', 'Yaquina Bay Monitoring Site');
geobasemap topographic
legend([a b c d], 'FontSize',16, 'Location', 'southeast')
grid off

