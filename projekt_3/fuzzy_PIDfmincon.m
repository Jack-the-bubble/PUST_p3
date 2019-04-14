% PID rozmyty

%ograniczenia nastaw PID
    function xPID = fuzzy_PIDfmincon(x0PID)
    global lreg;
    
    %x0PID - punkt startowy
    Kmin = 0;
    Kmax = 10;

    Timin = 0;
    Timax = 1000;

    Tdmin = 0;
    Tdmax = 1000;

    %funkcja celu (PID) -> Zad6PID (na razie)

    
    %ograniczenia nierownosciowe
    APID_part = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
    APID_zeros = zeros(6,3);
        
    APID = [APID_part APID_zeros APID_zeros];
    APID_row = APID; 
    
    bPID_part = [Kmax; Timax; Tdmax; -Kmin; -Timin; -Tdmin];
    bPID = bPID_part;
        
    for i = 2 : lreg
        APID_row = [APID_zeros APID_row(:, 1:end-3)];
        APID = [APID; APID_row];
        
        bPID = [bPID; bPID_part];
    end
    
%     %punkt startowy
%      x0PID = [1; 30; 5];

    %dodatkowe opcje
    optionsPID = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'notify');
    
    %dodatkowe opcje
    %optionsPID = optimoptions(@ga, 'Display', 'iter', 'MaxGenerations', 30);
    %optionsPID = optimoptions(@ga, 'Display', 'none', 'MaxGenerations', 10);
    
    %x = [x(1) x(2) x(3)] = [K Ti Td]
    xPID = fmincon(@p_Zad5PIDRozm, x0PID, APID, bPID, [], [], [], [], [], optionsPID);
    %xPID = ga(@p_Zad5PIDRozm, 9, APID, bPID, [], [], [], [], [], optionsPID);
end


