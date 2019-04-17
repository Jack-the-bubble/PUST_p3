%Projekt PUST
%Zadanie 6
%Funkcja obliczaj�ca b��d DMC

function error = p_Zad5DMCRozm(x)
    clear e u y U Y
    global lreg st D;
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
    
    clear U Y

    %REGULATOR DMC -----------------------------------------------------
    
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
    du = 0;

    ypast = 0.0; %poprzednia wartosc wyjscia
    upast = zeros(lreg,1); %poprzednia wartosc sterowania
    e = 0.0; %uchyb

    U = ones(lreg,iterNum)*Upp;
    Y = ones(iterNum,1)*Ypp;
    dUpast = zeros(lreg, D-1); %wektor przeszlych przyrostow sterowan


    %przypisanie odpowiedzi skokowej (znormalizowanej=
       
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
    
    for k = 7 : iterNum

        %ypast = y(k-1);

        Y(k) = symulacja_obiektu3y(Ukonc(k-5),Ukonc(k-6),Y(k-1),Y(k-2));
        y(k) = Y(k) - Ypp;

        e = yZad(k) - y(k);
        
        
        for i=1:lreg
            upast(i) = u(i,k-1);

%             ue(i) = ke(i)*e;
%             uu = Ku*dUpast;
%             du = ue-uu;

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
            
            if lreg == 2
                if i == 1
                    value = trapmf(Ukonc(k-1),[-1 -1 0 0.6]);
                end
                
                if i == 2
                    value = trapmf(Ukonc(k-1),[0 0.6 1 1]);
                end
            end
            
            if lreg == 3
                if i == 1
                    value = trapmf(Ukonc(k-1),[-1 -1 -0.2 0]);
                end
                
                if i == 2
                    value = trapmf(Ukonc(k-1),[-0.2 0 0.4 0.6]);
                end
                
                if i == 3
                    value = trapmf(Ukonc(k-1),[0.4 0.6 1 1]);
                end
            end
            

            
            Ukonc(k) = Ukonc(k) + U(i,k)*value;
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
    plot(Ukonc);
    
%         if lreg == 3
%         nazwa1 = sprintf('wykresy_txt/DMC_fuzzy/U__DMC_K=%g,%g,%g_Ti=%g,%g,%g_Td=%g,%g,%g_E=%g_.txt',K,Ti,Td,error);
%         nazwa2 = sprintf('wykresy_txt/DMC_fuzzy/Y__PID_K=%g,%g,%g_Ti=%g,%g,%g_Td=%g,%g,%g_E=%g_.txt',K,Ti,Td,error);
%     end
%     
%     if lreg == 2
%         nazwa1 = sprintf('wykresy_txt/PID_fuzzy/U__PID_K=%g,%g_Ti=%g,%g_Td=%g,%g_E=%g_.txt',K,Ti,Td,error);
%         nazwa2 = sprintf('wykresy_txt/PID_fuzzy/Y__PID_K=%g,%g_Ti=%g,%g_Td=%g,%g_E=%g_.txt',K,Ti,Td,error);
%     end
%         nazwa3 = 'wykresy_txt/PID_fuzzy/Yzad_single.txt';
% 
%  file = fopen(nazwa1, 'w');
%  A = [(1:iterNum);Ukonc'];
%  fprintf(file, '%4.3f %.3f \n',A);
%  fclose(file);
% 
%  file = fopen(nazwa2, 'w');
%  B = [(1:iterNum);Y'];
%  fprintf(file, '%4.3f %.3f \n',B);
%  fclose(file);
% 
%  file = fopen(nazwa3, 'w');
%  C = [(1:iterNum);(yZad+Ypp)'];
%  fprintf(file, '%4.3f %.3f \n',C);
%  fclose(file);
    end
