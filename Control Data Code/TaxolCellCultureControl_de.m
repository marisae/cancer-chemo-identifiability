function F = TaxolCellCultureControl_de(t,y,pars,theta)

% P proliferation rate
lam = pars(2);

% R to P transition
K   = pars(1)*100;
V0  = pars(3)*K;
L0  = (K^theta)/((V0^theta) + (K^theta));

Ncel = y(1)+y(2);
Lfac = ((K-Ncel)^theta)/((V0^theta) + ((K-Ncel)^theta));

aRP = 0.9*(1/L0);   % per day from Kim_PrlifQuies

% The differntial equations

F = [-lam*y(1) + aRP*y(2)*Lfac
    2*lam*y(1) - aRP*y(2)*Lfac];