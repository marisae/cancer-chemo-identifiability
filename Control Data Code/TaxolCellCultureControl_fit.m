function G = TaxolCellCultureControl_fit(theta,i)

% In vitro cell growth data for D-54Mg glioma cells
% from Terzis et al. Brit J Cancer 1997;75:1744
% From Bowman et al. Glia 1999;27:22
% glioma cell volume is 0.916 picoliters, 1 mm^3 = 1e6 pl or ~1.091 million
% cells

Time = [0      3      6      9     12     15    ]';        % days
Cell = [0.009  0.050  0.120  0.189  0.230  0.260]'*1091;   % thousands of cells
Cerr = [0.006  0.012  0.010  0.011  0.011  0.011]'*1091;   % thousands of cells

data = [Time  Cell];

% figure(1)  
% hold on
% set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
% errorbar(Time,Cell,Cerr,'rs','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12,'LineWidth',2)
% xlabel('Time, in days')
% ylabel('Number of cells, in thousands')

% dN/dt at t = 0
dNdt0 = (Cell(2)-Cell(1))/(Time(2)-Time(1));
N0    = Cell(1);

% lb for lambda = dNdt0/N0.

f = @(pars,tdata)TaxolCellCultureControl_soln(pars,tdata,dNdt0,N0,theta);

% Fitting

parsic = [3    1+dNdt0/N0    1];
lb     = [0    dNdt0/N0      0];
ub     = [Inf     Inf        Inf];

[parfit,resnorm] = lsqcurvefit(f,parsic,data(:,1),data(:,2),lb,ub);

% Reults of fit

K   = parfit(1)*100;
V0  = parfit(3)*K;
lam = parfit(2);

P0 = dNdt0/lam;
R0 = N0-P0;

error = sqrt(resnorm);

% plotting fit
Tfinal = Time(1):0.1:Time(end);
Ffinal = TaxolCellCultureControl_soln(parfit,Tfinal,dNdt0,N0,theta);

figure(1)
plot(Tfinal,Ffinal,'LineWidth',2)

Niter = 0:.01:0.26;
Niter2 = Niter*1091; 
Lfactor = ((K-Niter2).^theta)./((V0^theta) + ((K-Niter2).^theta));
L0  = (K^theta)/((V0^theta) + (K^theta));
% figure(i+1)
% hold on
% set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
% plot(Niter,Lfactor/L0,'LineWidth',2)

% Fit results
%---------------
G(1) = K;
G(2) = lam;
G(3) = V0;
G(4) = P0;
G(5) = R0;
G(6) = resnorm;
G(7) = parfit(3);






