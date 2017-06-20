clear all

% name    a0        ka        r0        d0        kd    
% parsic = [8.3170    8.0959    0.0582    1.3307  119.1363];

% Data from Terzis et al. Brit J Cancer 1997;75:1744.
% From Bowman et al. Glia 1999;27:22, glioma cell volume is 0.916
% picoliters, 1 mm^3 = 1e6 pl or ~1.091 million cells

% We fit the Tx data only to fit 6 parameters. 3 exponents previously
% determined by TaxolTreatment_fit.

Time = [0      3      6      9     12     15    ]';        % days

% Control data
Cell = [0.009  0.050  0.120  0.189  0.230  0.260]'*1091;   % thousands of cells
Cerr = [0.006  0.012  0.010  0.011  0.011  0.011]'*1091;   % thousands of cells

% 0.005 ug/ml Taxol
Cell005 = [0.009  0.047  0.089  0.149  0.198  0.219]'*1091;   % thousands of cells
Cerr005 = [0.006  0.013  0.010  0.011  0.013  0.010]'*1091;   % thousands of cells

% 0.010 ug/ml Taxol
Cell010 = [0.009  0.043  0.077  0.093  0.109  0.128]'*1091;   % thousands of cells
Cerr010 = [0.006  0.012  0.013  0.012  0.014  0.012]'*1091;   % thousands of cells

% 0.040 ug/ml Taxol
Cell040 = [0.009  0.025  0.047  0.054  0.076  0.085]'*1091;   % thousands of cells
Cerr040 = [0.005  0.010  0.010  0.011  0.010  0.010]'*1091;   % thousands of cells

% 0.100 ug/ml Taxol
Cell100 = [0.009  0.025  0.026  0.028  0.029  0.031]'*1091;   % thousands of cells
Cerr100 = [0.006  0.010  0.009  0.008  0.011  0.011]'*1091;   % thousands of cells

C005 = mean(Cell005);
C010 = mean(Cell010);
C040 = mean(Cell040);
C100 = mean(Cell100);

sigmasq = (mean([(Cerr005/C005); (Cerr010/C010); (Cerr040/C040); (Cerr100/C100)]))^2;
threshold = sigmasq*chi2inv(0.95,5);

% Reshape data
tdata = [Time;Time;Time;Time];
cdata = [Cell005/C005;Cell010/C010;Cell040/C040;Cell100/C100];

data = [tdata  cdata];

%%

% name    a0        ka        r0        d0        kd    
parsic = [          8.0959    0.0582    1.3307  119.1363];
lb     = 0.9*parsic;
ub     = 1.1*parsic;

iend   = 61;

for i = 1:1:iend
    
    alpha0(i) = (1 + 0.05*(i-1))*8.3170;
    
    f = @(pars,tdata)testa0_soln(pars,tdata,alpha0(i),C005,C010,C040,C100);
     
    % Fitting
    [parfit,resnorm(i)] = lsqcurvefit(f,parsic,data(:,1),data(:,2),lb,ub);
    
    kalpha(i) = parfit(1);
    rho0(i)   = parfit(2);
    delta0(i) = parfit(3);
    kdelta(i) = parfit(4);
    
    i
    
    parsic = parfit;
    lb     = 0.9*parsic;
    ub     = 1.1*parsic;

end

% name    a0        ka        r0        d0        kd    
parsic = [          8.0959    0.0582    1.3307  119.1363];
lb     = 0.9*parsic;
ub     = 1.1*parsic;

for i = iend+1:1:iend+14
    
    alpha0(i) = (1 - 0.05*(i-iend))*8.3170;
    
    f = @(pars,tdata)testa0_soln(pars,tdata,alpha0(i),C005,C010,C040,C100);
     
    % Fitting
    [parfit,resnorm(i)] = lsqcurvefit(f,parsic,data(:,1),data(:,2),lb,ub);
    
    kalpha(i) = parfit(1);
    rho0(i)   = parfit(2);
    delta0(i) = parfit(3);
    kdelta(i) = parfit(4);
    
    i
    
    parsic = parfit;
    lb     = 0.9*parsic;
    ub     = 1.1*parsic;

end

[blah In] = sort(alpha0);

%% Figures

profparam = 1;
fixparams = 1;
paramnames = {'alpha0','kalpha','rho0','delta0','kdelta'};
labelnames = {'\alpha_0','k_\alpha','\rho_0','\delta_0','k_{\delta}'};
otherparams = paramnames;
otherparams(fixparams) = [];
otherlabels = labelnames;
otherlabels(fixparams) = [];
i=1;



figure(1)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(alpha0(In)/alpha0(1),resnorm(In),'LineWidth',2)
plot(alpha0(In)/alpha0(1),(threshold+resnorm(1))*ones(length(alpha0),1),'LineWidth',2)
plot(alpha0(1)/alpha0(1),resnorm(1),'^','MarkerSize',12,'LineWidth',2)
xlabel('Fold change in \alpha_0, maximum rate of G2/M cell arrest')
ylabel('Residual error in best fit')
ylim([0 5])

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,strcat(paramnames{profparam},'_error'),'pdf')
saveas(fig,strcat(paramnames{profparam},'_error'),'fig')


figure(2)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(alpha0(In),kalpha(In),'s','MarkerSize',12,'LineWidth',2)
xlabel('\alpha_0, maximum rate of G2/M cell arrest')
ylabel('k_{\alpha}')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'pdf')
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'fig')
i=i+1;


figure(3)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(alpha0(In),rho0(In),'s','MarkerSize',12,'LineWidth',2)
xlabel('\alpha_0, maximum rate of G2/M cell arrest')
ylabel('\rho_0')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'pdf')
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'fig')
i=i+1;


figure(4)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(alpha0(In),delta0(In),'s','MarkerSize',12,'LineWidth',2)
xlabel('\alpha_0, maximum rate of G2/M cell arrest')
ylabel('\delta_0')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'pdf')
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'fig')
i=i+1;


figure(5)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(alpha0(In),kdelta(In),'s','MarkerSize',12,'LineWidth',2)
xlabel('\alpha_0, maximum rate of G2/M cell arrest')
ylabel('k_{\delta}')

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'pdf')
saveas(fig,strcat(paramnames{profparam},'_',otherparams{i}),'fig')
i=i+1;
