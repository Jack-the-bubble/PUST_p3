%% DMC
function xDMC = single_DMCga(
%ograniczenia parametr�w DMC
%Nmax = D;
Nmax = 182;
Nmin = 1;

%Numax = D;
Numin = 1;

lambdamax = 1000;
lambdamin = 0.1;

%funkcja celu (DMC) -> Zad6DMC (na razie)

%ograniczenia nierownosciowe
ADMC = [1 0 0; -1 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1];
bDMC = [Nmax; 0; lambdamax; -Nmin; -Numin; -lambdamin];

%punkt startowy
x0DMC = [60; 30; 5];

%dodatkowe opcje
optionsDMC = optimoptions(@ga, 'Display', 'iter', 'MaxGenerations', 100);

%N i Nu ca�kowite:
IntCon = [1 2];

%x = [x(1) x(2) x(3)] = [N Nu lambda]
%DMC = fmincon(@Zad6DMC, x0DMC, ADMC, bDMC);
xDMC = ga(@Zad6DMC, 3, ADMC, bDMC, [], [], [], [], [], IntCon, optionsDMC);