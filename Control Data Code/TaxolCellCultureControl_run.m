clear all

trange  = 1:1:15;
nt = length(trange);

% In vitro cell growth data for D-54Mg glioma cells
% from Terzis et al. Brit J Cancer 1997;75:1744
% From Bowman et al. Glia 1999;27:22
% glioma cell volume is 0.916 picoliters, 1 mm^3 = 1e6 pl or ~1.091 million
% cells

TimeTT = [0      3      6      9     12     15    ]';        % days
CellTT = [0.009  0.050  0.120  0.189  0.230  0.260]'*1091;   % thousands of cells
CerrTT = [0.006  0.012  0.010  0.011  0.011  0.011]'*1091;   % thousands of cells

sigmasq = (mean(CerrTT))^2;
threshold = sigmasq*chi2inv(0.95,3)

figure(1)  
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
errorbar(TimeTT,CellTT,CerrTT,'rs','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',12,'LineWidth',2)
xlabel('Time, in days')
ylabel('Number of cells, in thousands')


for i = 1:1:nt
    
    G = TaxolCellCultureControl_fit(trange(i),i);
    
    VolumK(i)   = G(1);
    Volumlam(i) = G(2);
    VolumV0(i)  = G(3);
    
    VolumP0(i)  = G(4);
    VolumR0(i)  = G(5);
    Volumerr(i) = G(6);
    
    Volum3(i) = G(7);
    
    clear G K lam P0 R0 error N0 dNdt0 F F1
    
end

figure(21)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(trange,VolumK,'LineWidth',2)
xlabel('\theta')
ylabel('Carrying capacity K')

figure(22)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(trange,Volumlam,'LineWidth',2)
xlabel('\theta')
ylabel('Proliferation rate \lambda')

figure(24)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(trange,VolumV0,'LineWidth',2)
xlabel('\theta')
ylabel('V_0')

figure(23)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(trange,Volumerr,'LineWidth',2)
plot(trange,(threshold+min(Volumerr))*ones(length(trange),1),'LineWidth',2)
xlabel('\theta')
ylabel('Error in fit')


[MinErr, I] = min(Volumerr)
theta = trange(I)
K = VolumK(I)
lam = Volumlam(I)
V0 = VolumV0(I)
P0 = VolumP0(I)
R0 = VolumR0(I)
%%
% plotting g1/s to g2/m transition
Ncells = 0:0.01:0.26;
Ncel = Ncells*1091;
L0   = (K^theta)/((V0^theta) + (K^theta));
Lfac = ((K-Ncel).^theta)./((V0^theta) + ((K-Ncel).^theta));
aRP  = 0.9*(1/L0);   % per day from Kim_PrlifQuies
StoM = aRP*Lfac;

figure(20)
hold on
set(gca,'LineWidth',1.25,'FontSize',24,'FontWeight','normal','FontName','Helvetica')
plot(Ncells,StoM,'LineWidth',2)
xlabel('Number of cells, in thousands')
ylabel('G1/S to G2/M transition rate, per day')


%%%% Minima at theta = 9
%%%% K = 890.2742
%%%% lam = 2.6787
%%%% V0 = 1.1683e+03
%%%% P0 = 5.5662
%%%% R0 = 4.2528
