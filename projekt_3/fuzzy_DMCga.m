% DMC rozmyty

%ograniczenia nastaw PID
    function xDMC = fuzzy_DMCga()
    global lreg;
    
    %ograniczenia parametr�w DMC
    %Nmax = D;
    Nmax = 53;
    Nmin = 1;

    %Numax = D;
    Numin = 1;

    lambdamax = 1000;
    lambdamin = 0.1;

    %ograniczenia nierownosciowe
    ADMC_part = [1 0 0; -1 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
    ADMC_zeros = zeros(6,3);
        
    ADMC = [ADMC_part ADMC_zeros ADMC_zeros];
    ADMC_row = ADMC; 
    
    bDMC_part = [Nmax; 0; lambdamax; -Nmin; -Numin; -lambdamin];
    bDMC = bDMC_part;
        
    for i = 2 : lreg
        ADMC_row = [ADMC_zeros ADMC_row(:, 1:end-3)];
        ADMC = [ADMC; ADMC_row];
        
        bDMC = [bDMC; bDMC_part];
    end
    
    %dodatkowe opcje
    %optionsPID = optimoptions(@ga, 'Display', 'iter', 'MaxGenerations', 30);
    optionsDMC = optimoptions(@ga, 'Display', 'none', 'MaxGenerations', 30);
    IntCon = [1 2 4 5 7 8];
    xDMC = ga(@p_Zad5DMCRozm, 9, ADMC, bDMC, [], [], [], [], [], IntCon, optionsDMC);
end


