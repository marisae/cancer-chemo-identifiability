clear all
% Data from Terzis et al. Brit J Cancer 1997;75:1744.
% From Bowman et al. Glia 1999;27:22, glioma cell volume is 0.916
% picoliters, 1 mm^3 = 1e6 pl or ~1.091 million cells

% We fit the Tx data only to fit 5 parameters. 3 exponents previously
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
%%
sigmasq = (mean([(Cerr005/C005); (Cerr010/C010); (Cerr040/C040); (Cerr100/C100)]))^2;
%%
% Dose response curve, (end points from above data)
S005 = Cell005(end)/Cell(end);
S010 = Cell010(end)/Cell(end);
S040 = Cell040(end)/Cell(end);
S100 = Cell100(end)/Cell(end);
Taxol= [0  0.005  0.010  0.040  0.100]';    %ug/ml of Taxol
Surv = [1  S005   S010   S040   S100 ]'; 

% Reshape data
tdata = [Time;Time;Time;Time];
cdata = [Cell005/C005;Cell010/C010;Cell040/C040;Cell100/C100];

data = [tdata  cdata];

figure(1)  
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
errorbar(Time,Cell,Cerr,'s','MarkerSize',12,'LineWidth',2)
errorbar(Time,Cell005,Cerr005,'s','MarkerSize',12,'LineWidth',2)
errorbar(Time,Cell010,Cerr010,'s','MarkerSize',12,'LineWidth',2)
errorbar(Time,Cell040,Cerr040,'s','MarkerSize',12,'LineWidth',2)
errorbar(Time,Cell100,Cerr100,'s','MarkerSize',12,'LineWidth',2)
xlabel('Time, in days')
ylabel('Number of cells, in thousands')
xlim([-0.1 15.1])

figure(2)  
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(Taxol,Surv,'s','MarkerSize',12,'LineWidth',2)
xlabel('Taxol dose, in ug/ml')
ylabel('Surviving fraction of cells')

% par     a0        ka        r0        d0        kd     
parsic = [8.3170    8.0959    0.0582    1.3307  119.1363];
lb     = 0.*parsic;
lb(3)  = 0.99*parsic(3);
ub     = 2*parsic;

% parsic = [8.3170    8.0959    0.0582    1.3307  119.1363];
% lb = 0.999*parsic;
% ub = 1.001*parsic;

f = @(pars,tdata)TaxolTx_soln(pars,tdata,C005,C010,C040,C100);

% Fitting
[parfit,resnorm, resid,eflag,OUTPUT,LAMBDA,jacob] = lsqcurvefit(f,parsic,data(:,1),data(:,2),lb,ub);

parfit
resnorm

% plotting fit
P0 = 7.2700;
R0 = 2.5490;
y0 = [P0  R0  0];

dose = [0  5  10  40  100];

Tfinal = Time(1):0.1:Time(end);

% Control
[t,G1] = ode23s(@TaxolTx_de,Tfinal,y0,[],parfit,dose(1));
    
figure(1)
plot(Tfinal,G1(:,1)+G1(:,2)+G1(:,3),'LineWidth',2)

SFctrl = (G1(end,1) + G1(end,2) + G1(end,3));
SF(1)  = SFctrl/SFctrl;

clear t G1

for i = 2:1:5
    
    [t,G1] = ode23s(@TaxolTx_de,Tfinal,y0,[],parfit,dose(i));
    
    figure(1)
    plot(Tfinal,G1(:,1)+G1(:,2)+G1(:,3),'LineWidth',2)
    
    SF(i)= (G1(end,1) + G1(end,2) + G1(end,3))/SFctrl;
    
    clear t G1
    
end

figure(1)
legend({'Experimental data, control','Experimental data, 0.005 \mug/ml taxol','Experimental data, 0.010 \mug/ml taxol','Experimental data, 0.040 \mug/ml taxol','Experimental data, 0.100 \mug/ml taxol','Model fit','Model fit','Model fit','Model fit','Model fit'},'FontSize',16,'Location','northwest')
legend('boxoff')

figure(2)
plot(Taxol,SF,'LineWidth',2)

figure(1)

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Taxol_Treat_Timecourse','pdf')
saveas(fig,'Taxol_Treat_Timecourse','fig')

figure(2)

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
saveas(fig,'Taxol_Treat_DoseResp','pdf')
saveas(fig,'Taxol_Treat_DoseResp','fig')


%Estimated parameter variance
%format shortE
jacob = full(jacob); 
fim   = (jacob'*jacob)/sigmasq;
rank(fim)
eig(fim)
(eig(fim))/(max(eig(fim)))
cond(fim)
varp = inv(fim);
stdp  = sqrt(diag(varp)); %standard deviation is square root of variance
stdp  = 100*stdp'./parfit %[%]
1./diag(fim)'






