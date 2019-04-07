%Projekt PUST
%Zadanie 6
%Funkcja obliczaj�ca b��d DMC

function error = p_Zad5DMCRozm(lreg, x, st)
    clear e u y U Y

    %dane
    iterNum = 1000;
    Umin = -1;
    Umax = 1;
    Ypp = 0;
    Upp = 0;
  
%     trajektoria zadana dla rozmytego
    yZad = zeros(iterNum, 1);
    chunk = 120;
    yZad(1:chunk)=-1*0.14;
    yZad(chunk + 1:2*chunk)=-1*0.08;
    yZad(2*chunk +1:3*chunk)=0.26;
    yZad(3*chunk + 1:4*chunk)=0.896;
    yZad(4*chunk+1:5*chunk)=3.576;
    yZad(5*chunk +1 :iterNum)=5.594;
    
    %trapezowe funkcje przynaleznosci
    mf = zeros(lreg,2001);
    ux = -1:0.001:1;
    mf(1,:) = trapmf(ux,[-1 -1 -0.2 0]); %-0.5
    mf(2,:) = trapmf(ux,[-0.2 0 0.4 0.6]); %0.2
    mf(3,:) = trapmf(ux,[0.4 0.6 1 1]); %0.7
    
    clear U Y

    %REGULATOR DMC -----------------------------------------------------
    %horyzonty
    D = 53;
    %N = 35;
    %Nu = 2;
    %lambda = 21;
    
    N = zeros(lreg,1);
    Nu = zeros(lreg,1);
    lambda = zeros(lreg,1);
    
    for i=1:lreg
        N(i) = x(i,1);
        Nu(i) = x(i,2);
        lambda(i) = x(i,3);
    end

    %PARAMETRY 
    y = zeros(iterNum,1);
    u = zeros(lreg,iterNum);
    du = 0;

    ypast = 0.0; %poprzednia wartosc wyjscia
    upast = zeros(lreg,1); %poprzednia wartosc sterowania
    e = 0.0; %uchyb

    U = ones(lreg,iterNum)*Upp;
    Y = ones(iterNum,1)*Ypp;
    dUpast = zeros(lreg, D-1); %wektor przeszlych przyrostow sterowan


    %przypisanie odpowiedzi skokowej (znormalizowanej=
       
    M = cell(lreg,1);
    Mp = cell(lreg,1);
    Ku = cell(lreg,1);
    ke = zeros(lreg,1);

    for v=1:lreg
        % Macierz M
        M{v,1}=zeros(N,Nu);
        for i=1:N
           for j=1:Nu
              if (i>=j)
                 M(i,j)=st(v,i-j+1);
              end
           end
        end

        % Macierz Mp
        Mp{v,1}=zeros(N,D-1);
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
        I=eye(Nu);
        K=((M'*M+lambda*I)^(-1))*M';
        Ku{v,1}=K(1,:)*Mp;
        ke(v,1)=sum(K(1,:));

    end
    % -------------- DO REGULACJI ---------------

    du = zeros(lreg,1);
    Ukonc = zeros(1,iterNum);
    
    for k = 7 : iterNum

        %ypast = y(k-1);

        Y(k) = symulacja_obiektu3y(U(k-5),U(k-6),Y(k-1),Y(k-2));
        y(k) = Y(k) - Ypp;

        e = yZad(k) - y(k);
        
        
        for i=1:lreg
            upast(i) = u(i,k-1);

%             ue(i) = ke(i)*e;
%             uu = Ku*dUpast;
%             du = ue-uu;

            du(i) = ke(i)*e - Ku{i,1}*dUpast(i);
            
            u(i,k) = upast(i)+du(i);

            U(k) = u(k) + Upp;

            if U(i,k) < Umin 
                 U(i,k) = Umin;
                 du(i) = Umin - U(i,k-1);

            elseif U(i,k) > Umax 
                 U(i,k) = Umax;
                 du(i) = Umax - U(i,k-1);

            end

            dUpast(i) = [du(i); dUpast(i,1:end-1)];
            
            %zaokraglenie Ukonc(k-1), aby odczytywac wartosci z funkcji
            %przynaleznosci
            Ukonc(k-1) = round(Ukonc(k-1),3);
            
            %mnozenie przez funkcje przynaleznosci
            %1000*Ukonc(k-1)+1001 to przerobienie wartosci sterowania na
            %indeksy w macierzy funkcji przynaleznosci
            Ukonc(k) = Ukonc(k) + U(i,k)*mf(i,1000*Ukonc(k-1)+1001);
        end
    end
    
    Ydmc = Y;
    Udmc = U;
    
    error = sum(((yZad+Ypp) - Ydmc).^2);
    disp('Wskaznik jakosci regulacji dla regulatora DMC: '+ error);
    
%     figure(1)
% plot(U); hold on; plot(Y); hold off;hold on; plot(yZad+Ypp); hold off;    
    figure(2)   
    subplot(2,1,1);
    plot(Y);
    hold on;
    plot(yZad+Ypp);
    hold off;
    title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
    legend('y','yzad')
    subplot(2,1,2);
    plot(U);
    end