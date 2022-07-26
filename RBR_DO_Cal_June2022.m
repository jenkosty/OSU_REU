% RBR_DO_Cal_June2022.m
% 29 June 2022

clear
more off

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% O2 Saturation Data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% O2 Saturation Data (temperature and corresponding concentration)
Tsat = [20:1:30] ;
DOsatmgL = [9.09 8.92 8.74 8.58 8.42 8.26 8.11 7.97 7.83 7.69 7.56] ;
DOsatmmL = (DOsatmgL/32)*1000 ; %%% Converting to umol/L

%%% Displaying O2 Saturation Data
figure();
hold on
plot(Tsat,DOsatmmL,'ko','MarkerFaceColor','k') ;
axis([19 31 235 285]) ;
xlabel('Temperature (^oC)');
ylabel('DO Saturation Concentration (mm/L)');

%%% Interpolating to fill temperature values
temp_interp = 20:0.1:30;
DO_interp = interp1(Tsat, DOsatmmL, temp_interp);
plot(temp_interp, DO_interp, 'r');

%%
%%%%%%%%%%%%%%%%%%%
%%% First Trial %%%
%%%%%%%%%%%%%%%%%%%

%%% load the data (first trial)
load 011721_20220627_1425.mat
T1 = RBR.data(:,2) ; %%% Temperature data
DO11 = RBR.data(:,1) ; %%% DO data
DO21 = RBR.data(:,12) ; %%% DO data (different unit)
tm1 = (1/6)*(0:1:length(T1)-1)' ; %%% Creating time vector

%%% Used to calculate e-folding time scale
nn = find((tm1>=750).*(tm1<=850)) ;
DO0 = mean(DO11(nn)) ;
t1 = 562 ;
t2 = 700 ;
nn = find((tm1>=t1).*(tm1<=t2)) ;

%%% Used to calculate scaling factor
t1_scaling = 435;
t2_scaling = 480;
mm = find((tm1>=t1_scaling).*(tm1<=t2_scaling));
mean_temp = mean(T1(mm));
mean_DO = mean(DO11(mm));

mean_DO_sat = interp1(Tsat, DOsatmmL, mean_temp);

disp('The Trial 1 scaling factor is ' + string(mean_DO_sat/mean_DO));

%%% Displaying data from Trial 1
figure('Renderer', 'painters', 'Position', [100 100 1000 700]);
sgtitle('First Trial')

%%% Temperature subplot
ax1 = subplot(2,1,1) ;
hold on
plot(tm1,T1,'k') ;
xline(t1, 'g')
xline(t2, 'r')
xline(t1_scaling, '--g')
xline(t2_scaling, '--r')
ylabel('Temperature (^oC)')

%%% DO subplot
ax2 = subplot(2,1,2) ;
hold on
plot(tm1,DO11,'k') ;
xline(t1, 'g')
xline(t2, 'r')
xline(t1_scaling, '--g')
xline(t2_scaling, '--r')
ylabel('DO (\mu mol/L)')
xlabel('time (s)')

linkaxes([ax1 ax2], 'x');

%%% Displaying data for e-folding time scale
figure();
plot(tm1(nn)-t1,log(DO11(nn)-DO0),'k')
hold on
x = tm1(nn)-t1;
%plot(x, log(DO11(nn(1)))-x./21, 'r')
p = polyfit(x, log(DO11(nn)-DO0), 1);
f = polyval(p, x);
plot(x, f, 'b');
xlabel('time (s)');

disp('The Trial 1 e-folding time scale is ' + string(-1/p(1)) + 's')


%%
%%%%%%%%%%%%%%%%%%%%
%%% Second Trial %%%
%%%%%%%%%%%%%%%%%%%%

%%% load the data (second trial)
load 011721_20220627_1447.mat
T1 = RBR.data(:,2) ;
DO11 = RBR.data(:,1) ;
DO21 = RBR.data(:,12) ;
tm1 = (1/6)*(0:1:length(T1)-1)' ;

%%% Used to find e-folding time scale
nn = find((tm1>=220).*(tm1<=320)) ;
DO0 = mean(DO11(nn)) ;
t1 = 92 ;
t2 = 250 ;
nn = find((tm1>=t1).*(tm1<=t2)) ;

%%% Used to calculate scaling factor
t1_scaling = 572;
t2_scaling = 613;
mm = find((tm1>=t1_scaling).*(tm1<=t2_scaling));
mean_temp = mean(T1(mm));
mean_DO = mean(DO11(mm));

mean_DO_sat = interp1(Tsat, DOsatmmL, mean_temp);

disp('The Trial 2 scaling factor is ' + string(mean_DO_sat/mean_DO));

%%% Displaying data from Trial 2
figure('Renderer', 'painters', 'Position', [100 100 1000 700]);
sgtitle('Second Trial')

%%% Temperature subplot
ax1 = subplot(2,1,1) ;
plot(tm1,T1,'k') ;
xline(t1, 'g')
xline(t2, 'r');
xline(t1_scaling, 'g--')
xline(t2_scaling, 'r--');
ylabel('Temperature (^oC)')

%%% DO subplot
ax2 = subplot(2,1,2) ;
hold on
plot(tm1,DO11,'k') ;
xline(t1, 'g')
xline(t2, 'r');
xline(t1_scaling, 'g--')
xline(t2_scaling, 'r--');
ylabel('DO (\mu mol/L)')
xlabel('time (s)')

linkaxes([ax1 ax2], 'x');

%%% Displaying data for e-folding time scale
figure();
plot(tm1(nn)-t1,log(DO11(nn)-DO0),'k')
hold on
x = tm1(nn)-t1;
%plot(x, log(DO11(nn(1)))-x./21, 'r')
p = polyfit(x, log(DO11(nn)-DO0), 1);
f = polyval(p, x);
plot(x, f, 'b');
xlabel('time (s)');

disp('The Trial 2 e-folding time scale is ' + string(-1/p(1)) + 's')
