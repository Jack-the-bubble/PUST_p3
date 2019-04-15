%PUST Lab3
%Regulatory DMC

FUZZY = 0;

if FUZZY == 0

    clear all
    lab1_zad3 = fullfile('Lab1Zad3b.m');

    run(lab1_zad3);

    %przypisanie odpowiedzi skokowej (znormalizowanej)
    st = Ynorm(2:length(Ynorm));
    
    figure
    plot(st);

    addpath('F:\SerialCommunication'); % add a path to the functions
    initSerialControl COM5 % initialise com port

    %dane
    iterNum = 1400;
    Umin = 0;
    Umax = 100;
    Ypp = 31.87;
    Upp = 28;

    %trajektoria zadana
    yZad = ones(iterNum, 1)*Ypp;
    yZad(351:700) = Ypp+5;
    yZad(701:1050) = Ypp+15;
    yZad = yZad - Ypp;

    %REGULATOR DMC -----------------------------------------------------
    %horyzonty
    D = 734;
    N = 150;
    Nu = 15;
    lambda = 0.1;

    %PARAMETRY 
    du = 0;
    upast = 0.0; %poprzednia wartosc sterowania
    e = 0.0; %uchyb

    u = zeros(iterNum, 1);
    U = ones(iterNum, 1)*Upp;
    Y = ones(iterNum, 1)*Ypp;
    dUpast = zeros(D-1, 1); %wektor przeszlych przyrostow sterowan

    % Macierz M
    M=zeros(N,Nu);
    for i=1:N
       for j=1:Nu
          if (i>=j)
             M(i,j)=st(i-j+1);
          end
       end
    end

    % Macierz Mp
    Mp=zeros(N,D-1);
    for i=1:N
       for j=1:D-1
          if (i+j)<=D-1
             Mp(i,j)=st(i+j)-st(j);
          else
             Mp(i,j)=st(D)-st(j);
          end      
       end
    end

    % Obliczanie parametrów regulatora
    I=eye(Nu);
    K=((M'*M+lambda*I)^(-1))*M';
    Ku=K(1,:)*Mp;
    ke=sum(K(1,:));

    % -------------- DO REGULACJI ---------------
    k = 2;
    while(1)
        upast = u(k-1);

        %pobranie wyjscia obiektu
        measurements = readMeasurements(1:7);
        Y(k)= measurements(1);
        y = Y(k)-Ypp;

        e = yZad(k) - y;

        ue = ke*e;
        uu = Ku*dUpast;

        du = ue-uu;
        u(k) = upast+du;
        U(k) = u(k)+Upp;

        if U(k) <  Umin 
             U(k) = Umin;
             du = Umin-U(k-1);

        elseif U(k) > Umax 
             U(k) = Umax;
             du = Umax - U(k-1);

        end

        dUpast = [du; dUpast(1:end-1)];   

        %sendControls([ 1, 2, 3, 4, 5, 6], ... send for these elements
        %                 [50 , 0, 0, 0, U(k), 0]);  % new corresponding control values
        
        sendNonlinearControls(U(k));
        
        %wykresy
        figure(1)
        subplot(2,1,1);
        plot(Y(1:k));
        hold on;
        plot(yZad(1:k)+Ypp);
        hold off;
        title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(U(1:k));
        drawnow;

        disp("U: " + U(k) +" Y: "+ Y(k) +" Yzad: "+ yZad(k)+Ypp);

        k = k+1;
        waitForNewIteration();
    end
    
else

    clear all;
    
    addpath('F:\SerialCommunication'); % add a path to the functions
    initSerialControl COM5 % initialise com port

    %dane
    iterNum = 1400;
    Umin = 0;
    Umax = 100;
    Ypp = 31.87;
    Upp = 28;

    %trajektoria zadana
    yZad = ones(iterNum, 1)*Ypp;
    yZad(351:700) = Ypp+5;
    yZad(701:1050) = Ypp+15;
    yZad = yZad - Ypp;
    
    clear U Y

    %REGULATOR DMC -----------------------------------------------------
    
    D = 
    lreg = 3;
    x = [D D 1; D D 1; D D 1];
    
    N = zeros(lreg,1);
    Nu = zeros(lreg,1);
    lambda = zeros(lreg,1);
    
    for i=1:lreg
        N(i) = x(3*(i-1)+1);
        Nu(i) = x(3*(i-1)+2);
        lambda(i) = x(3*(i-1)+3);
    end

    %PARAMETRY 
    y = zeros(iterNum,1);
    u = zeros(lreg,iterNum);

    ypast = 0.0; %poprzednia wartosc wyjscia
    upast = zeros(lreg,1); %poprzednia wartosc sterowania
    e = 0.0; %uchyb

    U = ones(lreg,iterNum)*Upp;
    Y = ones(iterNum,1)*Ypp;
    dUpast = zeros(lreg, D-1); %wektor przeszlych przyrostow sterowan

    %przypisanie odpowiedzi skokowej (znormalizowanej)
    %to jeszcze dojdzie
       
    %M = cell(lreg,1);
    %Mp = cell(lreg,1);
    Ku = cell(lreg,1);
    ke = zeros(lreg,1);

    for v=1:lreg
        % Macierz M
        M=zeros(N(v),Nu(v));
        for i=1:N
           for j=1:Nu
              if (i>=j)
                 M(i,j)=st(v,i-j+1);
              end
           end
        end

        % Macierz Mp
        Mp=zeros(N(v),D-1);
        for i=1:N
           for j=1:D-1
              if (i+j)<=D-1
                 Mp(i,j)=st(v,i+j)-st(v,j);
              else
                 Mp(i,j)=st(v,D)-st(v,j);
              end      
           end
        end

        % Obliczanie parametrów regulatora
        I=eye(Nu(v));
        K=((M'*M+lambda(v)*I)^(-1))*M';
        Ku{v,1}=K(1,:)*Mp;
        ke(v,1)=sum(K(1,:));

    end
    % -------------- DO REGULACJI ---------------

    du = zeros(lreg,1);
    Ukonc = zeros(1,iterNum);
    
    k = 2;
    while(1)

        %pobranie wyjscia obiektu
        measurements = readMeasurements(1:7);
        Y(k)= measurements(1);
        
        y(k) = Y(k) - Ypp;

        e = yZad(k) - y(k);
         
        for i=1:lreg
            upast(i) = u(i,k-1);

            du(i) = ke(i)*e - Ku{i,1}*dUpast(i,:)';
            
            u(i,k) = upast(i)+du(i);

            U(i,k) = u(i,k) + Upp;

            if U(i,k) < Umin 
                 U(i,k) = Umin;
                 du(i) = Umin - U(i,k-1);

            elseif U(i,k) > Umax 
                 U(i,k) = Umax;
                 du(i) = Umax - U(i,k-1);

            end

            dUpast(i,:) = [du(i), dUpast(i,1:end-1)];
            
            %zaokraglenie Ukonc(k-1), aby odczytywac wartosci z funkcji
            %przynaleznosci
            %Ukonc(k-1) = round(Ukonc(k-1),3);
            
            %mnozenie przez funkcje przynaleznosci
            %1000*Ukonc(k-1)+1001 to przerobienie wartosci sterowania na
            %indeksy w macierzy funkcji przynaleznosci
            %Ukonc(k) = Ukonc(k) + U(i,k)*mf(i,1000*Ukonc(k-1)+1001);
            
           

            if i == 1
                value = trapmf(Ukonc(k-1),[-1 -1 -0.2 0]);
            end

            if i == 2
                value = trapmf(Ukonc(k-1),[-0.2 0 0.4 0.6]);
            end

            if i == 3
                value = trapmf(Ukonc(k-1),[0.4 0.6 1 1]);
            end       

            Ukonc(k) = Ukonc(k) + U(i,k)*value;
        end
        
        %sendControls([ 1, 2, 3, 4, 5, 6], ... send for these elements
        %                 [50 , 0, 0, 0, U(k), 0]);  % new corresponding control values
        
        sendNonlinearControls(Ukonc(k));
        
        %wykresy
        figure(1)
        subplot(2,1,1);
        plot(Y(1:k));
        hold on;
        plot(yZad(1:k)+Ypp);
        hold off;
        title(['Regulator DMC N=',sprintf('%g',N'),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(Ukonc(1:k));
        drawnow;
        
        disp('U: ' + Ukonc(k) +' Y: '+ Y(k) +' Yzad: '+ yZad(k)+Ypp);
    
        k = k+1;
        waitForNewIteration();
        
    end
end