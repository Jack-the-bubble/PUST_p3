%PUST Projekt 1
%Zadanie 6
%Funkcja obliczaj�ca b��d PID

function error = p_Zad3PID(x)
    %nastawy i czas probkowania (TODO: pomy�le� o sekcjonowaniu kodu)
    %Ti = 24;
    %Td = 3;
    %K = 0.7;
    
    %nastawy
    K = x(1);
    Ti = x(2);
    Td = x(3);
    T = 0.5;
    
    %dane
    iterNum = 800;
    Umin = -1;
    Umax = 1;
    Ypp = 0;
    Upp = 0;
%     deltaUmax = 0.05;

%     yZad = ones(iterNum, 1)*2;
%     yZad = zeros(iterNum, 1);
%     chunk = 120
%     yZad(1:chunk)=3.576;
%     yZad(chunk + 1:2*chunk)=5.594;
%     yZad(2*chunk +1:3*chunk)=0.26;
%     yZad(3*chunk + 1:4*chunk)=0.896;
%     yZad(4*chunk+1:5*chunk)=-1*0.14;
%     yZad(5*chunk +1 :iterNum)=-1*0.08;



    yZad = zeros(iterNum, 1);
    chunk = 120;
    yZad(1:chunk)=-1*0.14;
    yZad(chunk + 1:2*chunk)=-1*0.08;
    yZad(2*chunk +1:3*chunk)=0.26;
    yZad(3*chunk + 1:4*chunk)=0.896;
    yZad(4*chunk+1:5*chunk)=3.576;
    yZad(5*chunk +1 :iterNum)=5.594;

    %wzorki
    r0 = K*(1+T/(2*Ti)+Td/T);
    r1 = K*(T/(2*Ti)-2*Td/T-1);
    r2 = K*Td/T;

    %sygna�y
    u = zeros(iterNum, 1);
    e = zeros(iterNum, 1);
    y = zeros(iterNum, 1);
    U = ones(iterNum, 1)*Upp;
    Y = ones(iterNum, 1)*Ypp;

    for k = 7 : iterNum
        %SYMULACJA ALGORYTMU
        %pobranie wyjscia obiektu
        Y(k)=symulacja_obiektu3y(U(k-5), U(k-6), Y(k-1), Y(k-2));

        %przesuniecie wyjscia i wart. zad. o punkt pracy
        y(k) = Y(k)-Ypp; 
    %     yZad(k) = yZad(k)-Ypp;

        %uchyb
        e(k) = yZad(k) - y(k);

        %sterowanie z przesunieciem o punkt pracy
        u(k) = r2*e(k-2)+r1*e(k-1)+r0*e(k)+u(k-1);
        U(k) = u(k)+Upp;

%         %ograniczenia na przyrosty sterowania
%         if U(k) - U(k-1) >= deltaUmax
%             U(k) = U(k-1) + deltaUmax;
%         elseif U(k) - U(k-1) <= -deltaUmax
%                 U(k) = U(k-1) - deltaUmax;
%         end

        %ograniczenia na wartosci sterowania
        if U(k) > Umax
            U(k) = Umax;
        elseif U(k) < Umin
            U(k) = Umin;
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
title(['Regulator PID K=',sprintf('%g',K'),' Ti=',sprintf('%g',Ti),' Td=',sprintf('%g',Td)]);
legend('y','yzad')
subplot(2,1,2);
plot(U);

 nazwa1 = sprintf('wykresy_txt/PID_single/U__PID_K=%g_Ti=%g_Td=%g_E=%g_.txt',K,Ti,Td,error);
 nazwa2 = sprintf('wykresy_txt/PID_single/Y__PID_K=%g_Ti=%g_Td=%g_E=%g_.txt',K,Ti,Td,error);
 nazwa3 = 'wykresy_txt/PID_single/Yzad_single.txt';

file = fopen(nazwa1, 'w');
A = [(1:iterNum);U'];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen(nazwa2, 'w');
B = [(1:iterNum);Y'];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

file = fopen(nazwa3, 'w');
C = [(1:iterNum);(yZad+Ypp)'];
fprintf(file, '%4.3f %.3f \n',C);
fclose(file);
end