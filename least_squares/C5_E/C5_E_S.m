%%% Sezione C5 est - SOLETTA (S)

%% Initialization
fprintf('INITIALIZATION...\n');
clc
clear all
close all
fclose all;
fontsize=14;
linethick=0.5;
format long

currentFolder = pwd;

%%
cs = 5/1000; % strain random noise amplification percentage on the mean value;
cT = 5/100; % temperature random noise amplification percentage on the mean value;

%% Format
timeformat='dd-mm-yyyy HH:MM:SS';
dateformat='dd-mm-yyyy';

% No. of days taken into account
maxNOofDays=1e5;

%% Importing data

fprintf('IMPORTING DATA...\n');

% % Windows
% Sdata=xlsread('..\data\FBG_C5_E_P.xlsx');
% Tdata = xlsread('..\data\T_C5_E_B0C.xlsx');


% Linux
Sdata=xlsread('../../data/FBG_C5_E_S.xlsx'); % Box girder
Tdata = xlsread('../../data/T_C5_E_B0C.xlsx'); % Column 5 (inner) & 6 (outer) for  thermocouple


%% Saving data to matlab format mat

fprintf('SAVING DATA TO mat FILE...\n');

save './S/Sdata.mat' Sdata;
save './S/Tdata.mat' Tdata;

%% Loading mat files

fprintf('LOADING mat FILES...\n');

load './S/Sdata.mat' Sdata;
load './S/Tdata.mat' Tdata;

%% Sorting data by time

fprintf('SORTING DATA BY TIME...\n');

Ssdata = sortrows(Sdata, 1);
Tsdata = sortrows(Tdata, 1);


%% Defining strain time and temperature time
% Time data
ts = Ssdata(:, 1);
tT = Tsdata(:,1);
%t=d(:,1)+693960; % 693960 to converto from Excel to MATLAB


%% Selecting same data time 

fprintf('SELECTING DATA WITH SAME TIME...\n');

[I1, I2] = find(abs(ts - tT.') <= 1/288); % time difference less than 5 mins
I1 = unique(I1);
I2 = unique(I2);

length(I1);
length(I2);

t1 = ts(I1);
t2 = tT(I2);

s = Ssdata(I1, 2);
Tinn = Tsdata(I2, 5);
Tout = Tsdata(I2, 6);

check = abs(t1 - t2) <= 1/288;

I = [];
n = 1;
for i = 1 : length(check)
    if check == 0
        I(n) = 0;
        n = n + 1;
    end
end
I

for i = 1 : length(Tout)
    if Tout(i) > 50
        Tout(i) = Tout(1);
    end
end

for i = 1 : length(Tinn)
    if Tinn(i) > 50
        Tinn(i) = Tinn(1);
    end
end


%% Purging NaN data
fprintf('PURGING NaN DATA...\n');


[tu, it] = find(t2 <= 737472);
t = t2(tu);
s = s(tu);
Tout = Tout(tu);
Tinn = Tinn(tu);



[tuu, itt] = find(t2 >= 736992 & t2 <= 737034);
Tout(tuu) = Tout(1);
Tinn(tuu) = Tinn(1);


T = Tinn;
% T = Tout;
% T = (Tout+Tinn)/2; % mean value

%% devstd strain gauges Colle Isarco

fprintf('CALCULATING RAW DATA MEAN STANDARD DEVIATION ON 5 DAYS...\n\n');

std1 = std(s(115:115+3));
std2 = std(s(115+24:115+24+3));
std3 = std(s(115+24*2:115+24*2+3));
std4 = std(s(115+24*3:115+24*3+3));
std5 = std(s(115+24*5:115+24*5+3));

stdmean = 1/5 * sqrt(std1^2+std2^2+std3^2+std4^2+std5^2);
fprintf('Raw data std mean: %f\n\n', stdmean);

save './S/export/stdmean.mat' stdmean;
writematrix(stdmean, './S/export/stdmean.csv')

%% Modified data
% Adding random error to data

fprintf('ADDING RANDOME NOISE TO DATA...\n');

% s = s + (rand(length(s), 1) - 1/2) .* s .* cs;
T = T + ([0; rand(length(T) - 1, 1)] - 1/2) .* T .* cT;

%% devstd strain gauges mod

fprintf('CALCULATING MODIFIED DATA MEAN STANDARD DEVIATION ON 5 DAYS...\n\n');

std1 = std(s(115:115+3));
std2 = std(s(115+24:115+24+3));
std3 = std(s(115+24*2:115+24*2+3));
std4 = std(s(115+24*3:115+24*3+3));
std5 = std(s(115+24*5:115+24*5+3));

stdmeanmod = 1/5 * sqrt(std1^2+std2^2+std3^2+std4^2+std5^2);
fprintf('Modified data std mean: %f\n\n', stdmeanmod);

save './S/export/stdmeanmod.mat' stdmeanmod;
writematrix(stdmeanmod, './S/export/stdmeanmod.csv')

%% Plot
fprintf('PLOTTING DATA...\n');

% Strain data

% figure('WindowState', 'normal')

figure(1);
set(gca,'FontSize',fontsize,'FontName','Times New Roman')

subplot(2,1,1)
hold on
plot(ts, Ssdata(:, 2) - Ssdata(1, 2),'-r','LineWidth',linethick);
xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\Delta\epsilon [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
datetick('x',dateformat);
title("Raw strain data");
grid on

subplot(2,1,2)
plot(t, s-s(1), '-b', 'LineWidth', linethick);
xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\Delta\epsilon [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
datetick('x',dateformat);
title("Selected strain data");
grid on
hold off

cd S/
% cd P/

% saveas(figure(1), 'raw/Fig_01_Strain_Comparison.jpg');

cd(currentFolder)

% Temperature data

figure('WindowState', 'normal')

figure(2);
set(gca,'FontSize',fontsize,'FontName','Times New Roman')

subplot(2,1,1)
hold on
plot(tT(7:end), Tsdata(7:end, 2) - Tsdata(7, 2),'-r','LineWidth',linethick);
xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\Delta T [°C]','FontSize',fontsize,'FontName','Times New Roman');
datetick('x',dateformat)
title("Raw temperature data");
grid on


subplot(2,1,2)
hold on
plot(t, Tout - Tout(1),'-b','LineWidth',linethick);
plot(t, Tinn - Tinn(1),'-r','LineWidth',linethick);
% plot(t, T - T(1),'-g','LineWidth',linethick);
legend("Thermocouple 5", "Thermocouple 4");
xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\Delta T [°C]','FontSize',fontsize,'FontName','Times New Roman');
datetick('x',dateformat)
title("Selected temperature data");
grid on
% 
% 
cd S/
% cd P/

% saveas(figure(2), 'raw/Fig_02_Temp_Comparison.jpg');

cd(currentFolder)


%% Compensation
% Daily compensation (season compensation needs alpha parameter)
    lastmeasure = t(end);
rec = datenum('01-04-2017 05:00:00', timeformat); % First expected measurement that has to be considered
n=1;
k=1;
I = [];
while ( rec < lastmeasure && k <= maxNOofDays )
    check = abs(t - rec); % Check returns the difference (in days) between each measurement and the one we are looking for
    [M, Itemp] = min(check); % M is the value of check. Itemp is the position of M in t array
    if M < (1/12) % The measurements considered have to be taken between 4 AM and 6 AM. 1/12 = 2 hrs
        I(n) = Itemp; % I is the vector of the indexes of the considered measurements
        n=n+1;
    end
    rec=rec+1; % increase of 1 day
    k=k+1;
end
n=n-1;
n
%% 
%n=50; %use this line if you want to set a fixed number of n

% ...purged are vectors containing the considered measurements, i.e. those taken around 5 AM
tpurged=t(I,:);
spurged=s(I,:);
Tpurged=T(I,:);

dataPurged = [tpurged, spurged, Tpurged];

%% Plot

fprintf('PLOTTING PURGED DATA...\n');

% Strain data

% figure('WindowState', 'normal')
figure(3);
set(gca,'FontSize',fontsize,'FontName','Times New Roman')

subplot(2,1,1)
hold on
plot(t, s - s(1),'-r','LineWidth',linethick);
xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\Delta\epsilon [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
datetick('x',dateformat)
title("Strain data");
grid on

subplot(2,1,2)
hold on
plot(tpurged, spurged - spurged(1),'-b','LineWidth',linethick);
xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\Delta\epsilon [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
datetick('x',dateformat)
title("Purged strain data");
grid on

cd("S/");
%
% saveas(figure(3), 'raw/Fig_03_Strain_Purged_Comparison.jpg');

cd(currentFolder)


%% Purge NaN data

fprintf('PURGING NaN DATA...\n');

k = 1;
for i = 1 : n
    if(isfinite(spurged(i)) && isfinite(Tpurged(i)))
            dt(k, 1) = tpurged(i) - tpurged(1);
            ds(k, 1) = spurged(i) - spurged(1);
            dT(k, 1) = Tpurged(i) - Tpurged(1);
            k = k+1;
    end
end
len = length(dt)


%% Parameter estimation
%Choosing number of days to analyze
fprintf('CHOOSING NUMBER OF DAYS...\n');

nd = len;
dtnd=dt(1:nd);
dsnd=ds(1:nd);
dTnd=dT(1:nd);
L=length(dtnd);

fprintf('number of days: %i\n', nd);

% LSA

fprintf('APPLYING LEAST SQUARE ANALYSIS...\n');

D = [ones(size(dtnd)), dtnd, dTnd];
ETHETA = D \ dsnd; %you do not need to compute the Pseudo Inverse Matrix. Matlab with the command \ calculates directly the regression coefficients. 
Y = dsnd - D * ETHETA; %residual
SIGMA = (Y' * Y) * (D' * D)^-1 / (length(dsnd) - length(ETHETA)); % covariation matrix
RO = corrcov(SIGMA); %correlation matrix

format short;

fprintf('PRINTIG RESULTS...\n\n');

fprintf('Valid measurements = %4.0f\n',L);
fprintf('E[epsilon0] = %3.5f;     E[m] = %3.5f;      E[alpha] = %3.5f\n',ETHETA(1,1),ETHETA(2,1),ETHETA(3,1)); % Expected values of the parameters
fprintf('std[epsilon0] = %3.5f;   std[m] = %3.5f;    std[alpha] = %3.5f\n\n',sqrt(SIGMA(1,1)),sqrt(SIGMA(2,2)),sqrt(SIGMA(3,3))); % Standard deviations of the parameters
RO

%% Saving parameters to file
fprintf('SAVING PARAMETERS TO FILE...\n');

% ETHETA
save './S/export/ETHETA.mat' ETHETA;
writematrix(ETHETA, './S/export/ETHETA.csv')

% COVARIATION MATRIX
save './S/export/SIGMA.mat' SIGMA;
writematrix(SIGMA, './S/export/SIGMA.csv')

% STD
STD = [sqrt(SIGMA(1,1)), sqrt(SIGMA(2,2)), sqrt(SIGMA(3,3))];
save './S/export/STD.mat' STD;
writematrix(STD, './S/export/STD.csv')

% CORRELATION MATRIX
save './S/export/RO.mat' RO;
writematrix(RO, './S/export/RO.csv')

%% Plotting compensated measurements

fprintf('PLOTTING COMPENSATED MEASUREMENTS...\n');

% figure('WindowState', 'normal')
figure(4);
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',fontsize,'FontName','Times New Roman')

subplot(2,1,1)
hold on
plot(tpurged(1:nd),spurged(1:nd)-spurged(1),'-b');
xlabel('t');
ylabel('\Delta\epsilon [\mu\epsilon]');
datetick('x',dateformat)
title('Purged strain data');
grid on       

subplot(2,1,2)
hold on
plot(tpurged(1:nd),(spurged(1:nd)-spurged(1))-ETHETA(3,1)*(Tpurged(1:nd)-Tpurged(1)),'-b');

plot(tpurged(1:nd),mymodel([tpurged(1:nd)-tpurged(1)  zeros(length(tpurged(1:nd)))],ETHETA),'-r'); % zeros instead temperature data
plot(tpurged(1:nd),mymodel([tpurged(1:nd)-tpurged(1)  zeros(length(tpurged(1:nd)))],ETHETA+[0; sqrt(SIGMA(2,2)); 0]),':r'); % error estimation above
plot(tpurged(1:nd),mymodel([tpurged(1:nd)-tpurged(1)  zeros(length(tpurged(1:nd)))],ETHETA-[0; sqrt(SIGMA(2,2)); 0]),':r'); % error estimation below

xlabel('t');
ylabel('\Delta\epsilon - E[\alpha]\DeltaT [\mu\epsilon]');
datetick('x',dateformat)
title('Purged strain data + linear regression');
grid on
% 

cd("S/");
% saveas(figure(4), 'raw/Fig_04_Strain_Purged_Regression_Comparison.jpg');
cd(currentFolder)
%
%%
figure(5)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',fontsize,'FontName','Times New Roman')

hold on
plot(tpurged(1:nd),(spurged(1:nd)-spurged(1))-ETHETA(3,1)*(Tpurged(1:nd)-Tpurged(1)),'-b');

plot(tpurged(1:nd),mymodel([tpurged(1:nd)-tpurged(1)  zeros(length(tpurged(1:nd)))],ETHETA),'-r'); % zeros instead temperature data
plot(tpurged(1:nd),mymodel([tpurged(1:nd)-tpurged(1)  zeros(length(tpurged(1:nd)))],ETHETA+[0; sqrt(SIGMA(2,2)); 0]),':r'); % error estimation above
plot(tpurged(1:nd),mymodel([tpurged(1:nd)-tpurged(1)  zeros(length(tpurged(1:nd)))],ETHETA-[0; sqrt(SIGMA(2,2)); 0]),':r'); % error estimation below

xlabel('t');
ylabel('\Delta\epsilon - E[\alpha]\DeltaT [\mu\epsilon]');
datetick('x',dateformat)
title('Purged strain data + linear regression');
grid on
%
cd("S/");
% saveas(figure(5), 'raw/Fig_05_Strain_Purged_Regression.jpg');
cd(currentFolder)