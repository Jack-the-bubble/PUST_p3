clear all;
global st;
global lreg;
global D;
D = 61;

DMC = 1;
FUZZY = 0;
lreg = 3;


if FUZZY == 0
    if DMC == 0
        %xPID = single_PIDfmincon([0.2000    4.0126    2.0669]);%E=104.0181
        %errorPID2 = p_Zad3PID(xPID);
        errorPID2 = p_Zad3PID([0.2000 4.0126 2.0669]);
        
    else
        %odp skokowa dla DMC na całej przestrzeni sterowania    
        czas_sym = 600;
        u = ones(czas_sym, 1);
        y=zeros(czas_sym, 1);

        for k = 7:czas_sym
            y(k) = symulacja_obiektu3y(u(k-5),u(k-6),y(k-1),y(k-2));
        end
%         plot(y);
%         xlabel('k');
%         ylabel('y(k)');
%         hold on;
%         normalizowanie odpowiedzi skokowej
        
        st = y(2:end);
        %errorDMC = p_Zad3DMC([53 53 1]);
        %nastawy: [53 46.0000 1.0000  119.5432] E=113.6145
        %xDMC = single_DMCga();
        errorDMC2 = p_Zad3DMC([46.0000 1.0000 119.5432]);
    end
    
else
    if DMC == 0
        if lreg == 2
            %xPID = fuzzy_PIDfmincon([0.4488 122.0959 1.1900 0.0892 3.3326 0.4033]);
            errorFPID = p_Zad5PIDRozm([0.4488  122.0959    1.1900    0.0892 3.3326    0.4033]);%błąd 93.5929 dla lreq=2
        end
        
        if lreg == 3
            %xPID = fuzzy_PIDfmincon([1.9083   93.9994    0.7405    0.3995 31.5524    0.0462    1.4088   77.6580 0.1715]);
            errorFPID2 =p_Zad5PIDRozm([4.3301   94.0455    0.6022    0.3664 29.5223    0.0372    1.5132   79.0692 0.1643]);%błąd 344.8967 dla lreq=3
        end

        
    else
        if lreg == 2
            st = p_odpskokFDMC();
           % xDMC = fuzzy_DMCga();
            errorFDMC = p_Zad5DMCRozm([61 61 100 61 61 100]);%błąd 158.9567
        end
        
        if lreg==3
            st = p_odpskokFDMC();
            %xDMC = fuzzy_DMCga();
            errorFDMC = p_Zad5DMCRozm([61 61 10 61 61 10 61 61 10]); % błąd 134.8527
        end
        
    end
end    
 






% nazwa1 = sprintf('sprawko_dane/DMC_bez_zak/U__DMC_D=%g_N=%g_Nu=%g_L=%g_E=%g_.txt',D,N,Nu,lambda,wskaznikDMC);
%     nazwa2 = sprintf('sprawko_dane/DMC_bez_zak/Y__DMC_D=%g_N=%g_Nu=%g_L=%g_E=%g_.txt',D,N,Nu,lambda,wskaznikDMC);
%     nazwa3 = 'sprawko_dane/DMC_bez_zak/Yzad.txt';

% file = fopen(nazwa1, 'w');
% A = [(1:czas_sym);u'];
% fprintf(file, '%4.3f %.3f \n',A);
% fclose(file);
% 
% file = fopen(nazwa2, 'w');
% B = [(1:czas_sym);y'];
% fprintf(file, '%4.3f %.3f \n',B);
% fclose(file);
% 
% file = fopen(nazwa3, 'w');
% C = [(1:czas_sym);yzad'];
% fprintf(file, '%4.3f %.3f \n',C);
% fclose(file);