function F1 = TaxolTx_soln(pars,tdata,C005,C010,C040,C100)

P0 = 7.2700;
R0 = 2.5490;

y0 = [P0  R0  0];

dose = [5  10  40  100];    % dose in ng/ml
datamean = [C005  C010  C040  C100];

Time = [0      3      6      9     12     15    ]';        % days

% treatment runs
for i = 1:1:4
    
    [t,F] = ode23s(@TaxolTx_de,Time,y0,[],pars,dose(i));
    G(:,i) = (F(:,1) + F(:,2) + F(:,3))/datamean(i);
    
    clear t F
    
end

F1 = reshape(G,24,1);