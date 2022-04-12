%% Initialization
fprintf('INITIALIZATION...\n');
clc
clear all
close all
fclose all;
fontsize=14;
linethick=0.5;
format long


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
Tdata = xlsread('../../data/T_C5_E_B0C.xlsx'); % Column 2 (inner) & 3 (outer) for  thermocouple


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
Tinn = Tsdata(I2, 2);
Tout = Tsdata(I2, 3);

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

%% Modified data
% Adding random error to data

% fprintf('ADDING RANDOME NOISE TO DATA...\n');
% 
% s = s + (rand(length(s), 1) - 1/2) .* s .* cs;
% T = T + [0; rand(length(T) - 1, 1)] .* abs(T) .* cT;

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

% saveas(figure(1), 'raw/Fig_01_Strain_Comparison.jpg');


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
% saveas(figure(2), 'raw/Fig_02_Temp_Comparison.jpg');


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

% saveas(figure(3), 'raw/Fig_03_Strain_Purged_Comparison.jpg');


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



%% Bayesian parameter estimation 
% Metropolis - Hastings logarithm

j=1;

X=[dT dt]; % X is a matrix containing Delta_Temperature in first column and delta_time in second column
Y=[ds]; %Y is a vector containing delta_strains
[ETHETA,SIGMA]= MH(X,Y, j); %ETHETA contains parameters estimation (e0, alpha, m, sigma_LH). SIGMA is covariation matrix
                                           % MH is Metropoli - Hastings
                                           % function.
                                           %Output = ETHETA, SIGMA
RO=corrcov(SIGMA); % correlation matrix

format short;

fprintf('  E[epsilon0] = %8.4f;     E[alpha] = %6.3f;     E[m] = %16.10f;       E[sigma_eps] = %6.3f\n', ETHETA(1,1),ETHETA(2,1),ETHETA(3,1),ETHETA(4,1)); % Expected values of the parameters
fprintf('std[epsilon0] = %8.4f;   std[alpha] = %6.3f;   std[m] = %16.10f;     std[sigma_eps] = %6.3f\n\n', sqrt(SIGMA(1,1)),sqrt(SIGMA(2,2)),sqrt(SIGMA(3,3)),sqrt(SIGMA(4,4))); % Standard deviations of the parameters
fprintf('Pearson''s correlation:')
RO
    
    % Plots of T-compensated measurements
    figure(300+j)
    set(gcf,'Color',[1 1 1])
    
    subplot(2,1,1)
    set(gca,'FontSize',fontsize,'FontName','Times New Roman')
    hold on
    plot(tpurged,spurged(:,j)-spurged(1,j),'-b','LineWidth',linethick);
    xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
    ylabel('\Delta\epsilon [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
    datetick('x',dateformat)
    grid on
    
    
    subplot(2,1,2)
    set(gca,'FontSize',fontsize,'FontName','Times New Roman')
    hold on
    plot(tpurged,spurged(:,j)-spurged(1,j)-ETHETA(2,1)*(Tpurged(:,j)-Tpurged(1,j)),'-b','LineWidth',linethick);
    plot(tpurged,mymodel([zeros(length(tpurged),1) tpurged-tpurged(1)],ETHETA(1:3,1)),'-r','LineWidth',linethick);
    plot(tpurged,mymodel([zeros(length(tpurged),1) tpurged-tpurged(1)],ETHETA(1:3,1)+[0; 0; sqrt(SIGMA(3,3))]),':r','LineWidth',linethick);
    plot(tpurged,mymodel([zeros(length(tpurged),1) tpurged-tpurged(1)],ETHETA(1:3,1)-[0; 0; sqrt(SIGMA(3,3))]),':r','LineWidth',linethick);
    xlabel('t','FontSize',fontsize,'FontName','Times New Roman');
    ylabel('\Delta\epsilon-E[\alpha]\DeltaT [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
    datetick('x',dateformat)
    grid on
    
%     saveas(figure(300+j),[num2str(300+j) '.fig']);
    
    
    % Plots prior and posterior distributions
    figure(400+j)
    set(gcf,'Color',[1 1 1])
    
    sigma=sqrt(SIGMA(1,1));
    mu=ETHETA(1,1);
    x=(mu-4*sigma):(8*sigma/100):(mu+4*sigma);
    subplot(2,2,1)
    set(gca,'FontSize',fontsize,'FontName','Times New Roman')
    hold on
    plot(x,normpdf(x,mu,sigma),'-b','LineWidth',linethick);
    plot(x,normpdf(x,0,20),'-k','LineWidth',linethick);
    xlabel('\epsilon_{0} [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
    ylabel('pdf','FontSize',fontsize,'FontName','Times New Roman');
    grid on
    
    sigma=sqrt(SIGMA(2,2));
    mu=ETHETA(2,1);
    x=(mu-4*sigma):(8*sigma/100):(mu+4*sigma);
    subplot(2,2,2)
    set(gca,'FontSize',fontsize,'FontName','Times New Roman')
    hold on
    plot(x,normpdf(x,mu,sigma),'-b','LineWidth',linethick);
    plot(x,normpdf(x,12,4),'-k','LineWidth',linethick);
    xlabel('\alpha [\mu\epsilon/°C]','FontSize',fontsize,'FontName','Times New Roman');
    ylabel('pdf','FontSize',fontsize,'FontName','Times New Roman');
    grid on
    
    sigma=sqrt(SIGMA(3,3));
    mu=ETHETA(3,1);
    x=(mu-4*sigma):(8*sigma/100):(mu+4*sigma);
    subplot(2,2,3)
    set(gca,'FontSize',fontsize,'FontName','Times New Roman')
    hold on
    plot(x,normpdf(x,mu,sigma),'-b','LineWidth',linethick);
    plot(x,normpdf(x,0,10/365),'-k','LineWidth',linethick);
    xlabel('m [\mu\epsilon/d]','FontSize',fontsize,'FontName','Times New Roman');
    ylabel('pdf','FontSize',fontsize,'FontName','Times New Roman');
    grid on
    
    sigma=sqrt(SIGMA(4,4));
    mu=ETHETA(4,1);
    x=(mu-4*sigma):(8*sigma/100):(mu+4*sigma);
    subplot(2,2,4)
    set(gca,'FontSize',fontsize,'FontName','Times New Roman')
    hold on
    plot(x,normpdf(x,mu,sigma),'-b','LineWidth',linethick);
    plot(x,normpdf(x,10,5),'-k','LineWidth',linethick);
    xlabel('\sigma_{LH} [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
    ylabel('pdf','FontSize',fontsize,'FontName','Times New Roman');
    grid on
    
%     saveas(figure(400+j),[num2str(400+j) '.fig']);