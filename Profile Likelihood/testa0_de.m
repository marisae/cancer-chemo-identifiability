function F = testa0_de(t,y,pars,drug,a0)

ka = pars(1);

r0 = pars(2);

d0 = pars(3);
kd = pars(4);

% Following parameters fit by running TaxolCellCulture_run
K   = 10.515*100;
V0  = 1.3907*K;
lam = 9.5722;

theta = 10;

% Values taken from 
aRP  = 20;     % per day from Kim_PrlifQuies

Ncel = y(1)+y(2)+y(3);
Lfac = ((K-Ncel)^theta)/((V0^theta) + ((K-Ncel)^theta));

% if running TaxolTreatment_fit, these exponenets are pars(7), pars(8) and
% pars(9), and each fit in turn. 
arstexp = 3;
adthexp = 4;

arst = a0*(drug^arstexp)/(ka^arstexp + (drug^arstexp));
adth = d0*(drug^adthexp)/(kd^adthexp + (drug^adthexp));
arcv = r0;

% The differntial equations

F = [-lam*y(1) + aRP*y(2)*Lfac - arst*y(1) + arcv*y(3)
     2*lam*y(1) - aRP*y(2)*Lfac
     arst*y(1) - adth*y(3) - arcv*y(3)];
 
 
 


