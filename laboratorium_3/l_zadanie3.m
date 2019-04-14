%Rozmyty PID

%PUST Lab3
%Zadanie 3

FUZZY = 0;

addpath('F:\SerialCommunication'); % add a path to the functions
initSerialControl COM5 % initialise com port

if (FUZZY == 0)
    %Nierozmyty PID
    
    %nastawy i czas probkowania
    K = 10;
    Ti = 60;
    Td = 0.25;
    
    T = 1;
    
    %dane
    iterNum = 1400;
    Umin = 0;
    Umax = 100;
    Ypp = 31.87;
    Upp = 28;

    yZad = ones(iterNum, 1)*Ypp;
    yZad(351:700) = Ypp+5; %potem zmienic na 50
    yZad(701:1050) = Ypp+15;
    yZad = yZad - Ypp;

    %wzorki
    r0 = K*(1+T/(2*Ti)+Td/T);
    r1 = K*(T/(2*Ti)-2*Td/T-1);
    r2 = K*Td/T;

    %sygnaly
    u = zeros(iterNum, 1);
    e = zeros(iterNum, 1);
    U = ones(iterNum, 1)*Upp; %zmienic tak, zeby bral od poczatku aktualne wartosci u i y
    Y = ones(iterNum, 1)*Ypp;
    
    k = 3;

    while(1)
        %SYMULACJA ALGORYTMU
        %pobranie wyjscia obiektu
        measurements = readMeasurements(1:7);
        Y(k)= measurements(1);

        %przesuniecie wyjscia i wart. zad. o punkt pracy
        y = Y(k)-Ypp; 
    %     yZad(k) = yZad(k)-Ypp;

        %uchyb
        e(k) = yZad(k) - y;

        %sterowanie z przesunieciem o punkt pracy
        u(k) = r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
        U(k) = u(k)+Upp;

        %ograniczenia na wartosci sterowania
        if U(k) > Umax
            U(k) = Umax;
        elseif U(k) < Umin
            U(k) = Umin;
        end
        
        %sendControls([ 1, 2, 3, 4, 5, 6], ... send for these elements
        %             [50 , 0, 0, 0, U(k), 0]);  % new corresponding control values
        sendNonlinearControls(U(k));
        
        
        %wykresy
        figure(1)
        subplot(2,1,1);
        plot(Y(1:k));
        hold on;
        plot(yZad(1:k)+Ypp);
        hold off;
        title(['Regulator PID K=',sprintf('%g',K'),' Ti=',sprintf('%g',Ti),' Td=',sprintf('%g',Td)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(U(1:k));
        drawnow;
        
        disp('U: ' + U(k) +' Y: '+ Y(k) +' Yzad: '+ yZad(k)+Ypp);
    
        k = k+1;
        waitForNewIteration();
        
    
    end
else
    % Rozmyty PID
    
    K = [10 10 10];
    Ti = [60 60 60];
    Td = [0.25 0.25 0.25];
    r0 = zeros(3,1);
    r1 = zeros(3,1);
    r2 = zeros(3,1);
    
    T = 1;
    
    for i=1:3
        r0(i) = K(i)*(1+T/(2*Ti(i))+Td(i)/T);
        r1(i) = K(i)*(T/(2*Ti(i))-2*Td(i)/T-1);
        r2(i) = K(i)*Td(i)/T;        
    end
    
    %dane
    iterNum = 1400;
    Umin = 0;
    Umax = 100;
    Ypp = 31.87;
    Upp = 28;

    yZad = ones(iterNum, 1)*Ypp;
    yZad(351:700) = Ypp+5;
    yZad(701:1050) = Ypp+15;
    yZad = yZad - Ypp;
    
    %trapezowe funkcje przynaleznosci
    mf = zeros(3,10001);
    ux = 0:0.01:100;
    %mf(1,:) = trapmf(ux,[-1 -1 -0.2 0]);
    %mf(2,:) = trapmf(ux,[-0.2 0 0.4 0.6]);
    %mf(3,:) = trapmf(ux,[0.4 0.6 1 1]);
    
    % Do zestawienia funkcji przynaleznosci z charakterystyka statyczna
%     figure(2)
%     hold on
%     plot(ux,mf);
    
    %sygnaly
    u = zeros(3, iterNum);
    e = zeros(iterNum, 1);
    y = zeros(iterNum, 1);
    U = ones(3, iterNum)*Upp;
    Y = ones(iterNum, 1)*Ypp;
    
    %ostateczne sterowanie po polaczeniu wszystkich regulatorow
    Ukonc = zeros(iterNum,1);

    k = 3;
    %SYMULACJA ALGORYTMU
    while(1)
        %pobranie wyjscia obiektu
        measurements = readMeasurements(1:7);
        Y(k)= measurements(1);
        
        %przesuniecie wyjscia i wart. zad. o punkt pracy
        y(k) = Y(k)-Ypp; 
    %     yZad(k) = yZad(k)-Ypp;

        %uchyb
        e(k) = yZad(k) - y(k);

        %sterowanie z przesunieciem o punkt pracy
        for j=1:3
            u(j,k) = r2(j)*e(k-2)+r1(j)*e(k-1)+r0(j)*e(k)+u(j,k-1);
            U(j,k) = u(j,k)+Upp;
        
            %ograniczenia na wartosci sterowania
            if U(j,k) > Umax
                U(j,k) = Umax;
            elseif U(j,k) < Umin
                U(j,k) = Umin;
            end
            
            %zaokraglenie Ukonc(k-1), aby odczytywac wartosci z funkcji
            %przynaleznosci
            %Ukonc(k-1) = round(Ukonc(k-1),2);
            %mnozenie przez funkcje przynaleznosci
            %20*Ukonc(k-1)+1 to przerobienie wartosci sterowania na
            %indeksy w macierzy funkcji przynaleznosci
            %Ukonc(k) = Ukonc(k) + U(j,k)*mf(j,100*Ukonc(k-1)+1);
            
            if i == j
                value = trapmf(Ukonc(k-1),[0 0 40 45]);
            end
            
            if i == j
                value = trapmf(Ukonc(k-1),[40 45 55 60]);
            end
            
            if i == j
                value = trapmf(Ukonc(k-1),[55 60 100 100]);
            end
            
            Ukonc(k) = Ukonc(k) + U(j,k)*value;
            
        end
        
        %sendControls([ 1, 2, 3, 4, 5, 6], ... send for these elements
        %             [50 , 0, 0, 0, U(k), 0]);  % new corresponding control values
        sendNonlinearControls(Ukonc(k));
        
        %wykresy
        figure(1)
        subplot(2,1,1);
        plot(Y(1:k));
        hold on;
        plot(yZad(1:k)+Ypp);
        hold off;
        title(['Regulator PID K=',sprintf('%g',K'),' Ti=',sprintf('%g',Ti),' Td=',sprintf('%g',Td)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(U(1:k));
        drawnow;
        
        disp('U: ' + Ukonc(k) +' Y: '+ Y(k) +' Yzad: '+ yZad(k)+Ypp);
    
        k = k+1;
        waitForNewIteration();

    end
end