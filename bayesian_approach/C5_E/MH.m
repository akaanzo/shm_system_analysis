function [ETHETA,SIGMA]=MH(X,Y,j)

format short;

global currentFolder;

N=10000; % number of samples
N1=N/10; 

theta0=[0;0;0;0]; % initial parameters value (row vector: {e0; alpha; m; sigma_LH})
% sigma=diag([2,1,10/365,1].^2) % error diagonal matrix (sigma_e0 = 2, sigma_alpha = 1, sigma_m = 10/365, sigma_LH = 1)
% sigma=diag([2,0.12704,0.00598,1.529527].^2) 
sigma=diag([2,0.12,3/365,1].^2); 
% sigma=diag([1,1,1,1].^2)

sample     =zeros(length(theta0),N); % matrix of random sample. Foreach parameter N random sample need to be created. Row_num = param_num. Col_num = N
sample(:,1)=theta0; % program initialization: first column is theta0 parameters
Pr         =logposterior(sample(:,1),X,Y); % posterior log-distribution using "logposterior" function created. First column will be created
i=2;

% AR logarithm
while i <= N
    candidate=mvnrnd(sample(:,i-1)',sigma)'; % mvnrnd generate random number gaussian distribution (candidate = theta'). Is a 4 elements vector
      if (candidate(4) > 0) % std dev must be greater than zero (sigma_LH > 0)
        PrCand=logposterior(candidate,X,Y); % candidate log_probability using logposterior(candidate, X, Y) calculated on candidate matrix [f(theta')]
        ap=min([1, (exp(PrCand-Pr))]); % accept probability = min{1, exp{log(PrCand/Pr)}} = min{1, exp{PrCand - Pr}}
      else
          ap=0; % if sigma_LH <= 0 candidate is rejected
      end
    u=rand; % randomize u in [0, 1]
    if (ap > u) %accept
        sample(:,i)=candidate; % sample ith column is accepted
        Pr=PrCand;
    else %reject
        sample(:,i)=sample(:,i-1); % keep previous sample ((i-1)th sample)
    end
    i=i+1;
end

samples=sample(:,N1:N); % sample matrix is matrix with column from N1 to N (cut first ~ 10%)
ETHETA=mean(samples')'; % posterior mean (mean of samples rows)
SIGMA=cov(samples'); % std dev

% ETHETA=mean(sample')';
% SIGMA=cov(sample');

%%

fontsize=14;
linethick=0.5;

figure(100+j)
set(gcf,'Color',[1 1 1])

subplot(4,1,1)
set(gca,'FontSize',fontsize,'FontName','Times New Roman')
hold on
plot(1:size(sample,2),sample(1,:),'-b','LineWidth',linethick);
xlabel('q','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\epsilon_{0} [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
xlim([0 N]);
grid on

subplot(4,1,2)
set(gca,'FontSize',fontsize,'FontName','Times New Roman')
hold on
plot(1:size(sample,2),sample(2,:),'-b','LineWidth',linethick);
xlabel('q','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\alpha [\mu\epsilon/°C]','FontSize',fontsize,'FontName','Times New Roman');
xlim([0 N]);
grid on

subplot(4,1,3)
set(gca,'FontSize',fontsize,'FontName','Times New Roman')
hold on
plot(1:size(sample,2),sample(3,:),'-b','LineWidth',linethick);
xlabel('q','FontSize',fontsize,'FontName','Times New Roman');
ylabel('m [\mu\epsilon/d]','FontSize',fontsize,'FontName','Times New Roman');
xlim([0 N]);
grid on

subplot(4,1,4)
set(gca,'FontSize',fontsize,'FontName','Times New Roman')
hold on
plot(1:size(sample,2),sample(4,:),'-b','LineWidth',linethick);
xlabel('q','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\sigma_{LH} [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
xlim([0 N]);
grid on

cd(currentFolder);
cd("S/raw/");

saveas(figure(100+j),[num2str(100+j) '.jpg']);

figure(200+j)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',fontsize,'FontName','Times New Roman')
plot3(sample(1,:),sample(2,:),sample(3,:),'-+','MarkerEdgeColor','b','MarkerSize',linethick*6);
xlabel('\epsilon0 [\mu\epsilon]','FontSize',fontsize,'FontName','Times New Roman');
ylabel('\alpha [\mu\epsilon/°C]','FontSize',fontsize,'FontName','Times New Roman');
zlabel('m [\mu\epsilon/d]','FontSize',fontsize,'FontName','Times New Roman');
grid on
saveas(figure(200+j),[num2str(200+j) '.jpg']);

cd(currentFolder);

end

