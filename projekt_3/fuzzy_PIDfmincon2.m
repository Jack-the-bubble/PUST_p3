% PID rozmyty

%ograniczenia nastaw PID
    function xPID = fuzzy_PIDfmincon2(x0PID)
    global lreg;
    
    %x0PID - punkt startowy
    Kmin = 0.1;
    Kmax = 10;

    Timin = 0.1;
    Timax = 1000000;

    Tdmin = 0;
    Tdmax = 10000;

    %funkcja celu (PID) -> Zad6PID (na razie)

    
    %ograniczenia nierownosciowe
    APID_part = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
    
    bPID_part = [Kmax; Timax; Tdmax; -Kmin; -Timin; -Tdmin];
    bPID = bPID_part;
    APID = APID_part
    for i = 1: lreg-1
        APID = [APID APID_part]
%         bPID = [bPID;  bPID_part]
    end
    
%     %punkt startowy
      x0PID = [0.5; 100; 0; 0.5; 100; 0; 0.5; 100; 0];

    %dodatkowe opcje
    optionsPID = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'iter');

    %x = [x(1) x(2) x(3)] = [K Ti Td]
    xPID = fmincon(@p_Zad5PIDRozm2, x0PID, APID, bPID, [], [], [], [], [], optionsPID);
    print('haha');
end


