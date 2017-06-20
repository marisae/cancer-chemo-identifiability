function F1 = TaxolCellCultureControl_soln(pars,tdata,dNdt0,N0,theta)

lam = pars(2);

P0 = dNdt0/lam;
R0 = N0-P0;

y0 = [P0  R0];

[t,F] = ode23s(@TaxolCellCultureControl_de,tdata,y0,[],pars,theta);

F1(:,1) = (F(:,1) + F(:,2));