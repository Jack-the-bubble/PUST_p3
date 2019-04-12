%% PID nierozmyty

%ograniczenia nastaw PID
    function xPID = single_PIDfmincon(x0PID)
    %x0PID - punkt startowy
    Kmin = 0.1;
    Kmax = 10;

    Timin = 0.1;
    Timax = 1000000;

    Tdmin = 0;
    Tdmax = 10000;

    %funkcja celu (PID) -> Zad6PID (na razie)

    %ograniczenia nierownosciowe
    APID = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
    bPID = [Kmax; Timax; Tdmax; -Kmin; -Timin; -Tdmin];

%     %punkt startowy
%      x0PID = [1; 30; 5];

    %dodatkowe opcje
    optionsPID = optimoptions(@fmincon, 'Algorithm', 'sqp', 'Display', 'iter');

    %x = [x(1) x(2) x(3)] = [K Ti Td]
    xPID = fmincon(@p_Zad3PID, x0PID, APID, bPID, [], [], [], [], [], optionsPID);
end


