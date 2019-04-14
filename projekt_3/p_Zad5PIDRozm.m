%PUST Projekt 1
%Zadanie 6
%Funkcja obliczaj�ca b��d PID

function error = p_Zad5PIDRozm(x)
    %nastawy i czas probkowania (TODO: pomy�le� o sekcjonowaniu kodu)
    %Ti = 24;
    %Td = 3;
    %K = 0.7;
    
    global lreg;
    
    K = zeros(lreg,1);
    Ti = zeros(lreg,1);
    Td = zeros(lreg,1);
    r0 = zeros(lreg,1);
    r1 = zeros(lreg,1);
    r2 = zeros(lreg,1);
    
    %nastawy
    T = 0.5;
    for i=1:lreg
        K(i) = x(3*(i-1)+1);
        Ti(i) = x(3*(i-1)+2);
        Td(i) = x(3*(i-1)+3);
        r0(i) = K(i)*(1+T/(2*Ti(i))+Td(i)/T);
        r1(i) = K(i)*(T/(2*Ti(i))-2*Td(i)/T-1);
        r2(i) = K(i)*Td(i)/T;        
    end
    
    %dane
    iterNum = 1270;
    Umin = -1;
    Umax = 1;
    Ypp = 0;
    Upp = 0;

    yZad = zeros(iterNum, 1);
    chunk = 120;
    yZad(1:chunk)=-1*0.14;
    yZad(chunk + 1:2*chunk)=-1*0.08;
    yZad(2*chunk +1:3*chunk)=0.26;
    yZad(3*chunk + 1:4*chunk)=0.896;
    yZad(4*chunk+1:5*chunk)=3.576;
    yZad(5*chunk +1 :iterNum)=5.594;
    
    %trapezowe funkcje przynaleznosci
%     mf = zeros(lreg,2001);
%     ux = -1:0.001:1;
%     mf(1,:) = trapmf(ux,[-1 -1 -0.2 0]);
%     mf(2,:) = trapmf(ux,[-0.2 0 0.4 0.6]);
%     mf(3,:) = trapmf(ux,[0.4 0.6 1 1]);
    
    % Do zestawienia funkcji przynaleznosci z charakterystyka statyczna
%     figure(2)
%     hold on
%     plot(ux,mf);
    
    %sygnaly
    u = zeros(lreg, iterNum);
    e = zeros(iterNum, 1);
    y = zeros(iterNum, 1);
    U = ones(lreg, iterNum)*Upp;
    Y = ones(iterNum, 1)*Ypp;
    
    %ostateczne sterowanie po polaczeniu wszystkich regulatorow
    Ukonc = zeros(iterNum,1);

    %SYMULACJA ALGORYTMU
    for k = 7 : iterNum
        %pobranie wyjscia obiektu
        Y(k)=symulacja_obiektu3y(Ukonc(k-5), Ukonc(k-6), Y(k-1), Y(k-2));

        %przesuniecie wyjscia i wart. zad. o punkt pracy
        y(k) = Y(k)-Ypp; 
    %     yZad(k) = yZad(k)-Ypp;

        %uchyb
        e(k) = yZad(k) - y(k);

        %sterowanie z przesunieciem o punkt pracy
        for j=1:lreg
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
            %Ukonc(k-1) = round(Ukonc(k-1),3);
            
            %mnozenie przez funkcje przynaleznosci
            %1000*Ukonc(k-1)+1001 to przerobienie wartosci sterowania na
            %indeksy w macierzy funkcji przynaleznosci
            %Ukonc(k) = Ukonc(k) + U(j,k)*mf(j,1000*Ukonc(k-1)+1001);
            %w = trapezoid(j,Ukonc(k-1));
            
            if j == 1
                value = trapmf(U(k-1),[-1 -1 -0.2 0]);
            end
    
            if j == 2
                value = trapmf(U(k-1),[-0.2 0 0.4 0.6]);
            end
    
            if j == 3
                value = trapmf(U(k-1),[0.4 0.6 1 1]);
            end
            
            Ukonc(k) = Ukonc(k) + U(j,k)*value;
        end

    end

    Ypid = Y;
    Upid = U;

     error = sum(((yZad+Ypp) - Ypid).^2);
     disp('Wskaznik jakosci regulacji dla regulatora PID: '+ error);
    
    figure(1)
    subplot(2,1,1);
    plot(Y);
    hold on;
    plot(yZad+Ypp);
    hold off;
    legend('y','yzad')
    subplot(2,1,2);
    plot(Ukonc);
end


function value = trapezoid(reg, u)
    if reg == 1
        value = trapmf(u,[-1 -1 -0.2 0])
    end
    
    if reg == 2
        value = trapmf(u,[-0.2 0 0.4 0.6])
    end
    
    if reg == 3
        value = trapmf(u,[0.4 0.6 1 1])
    end
end
