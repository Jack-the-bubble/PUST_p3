clear all;
DMC = 0;
FUZZY = 1;

global st;
global lreg;

lreg = 3;

if FUZZY == 0
    if DMC == 0
%          errorPID = p_Zad3PID([0.5 100 0]);
        %nastawy: [0.5 100 0] E=700,2026
        xPID = single_PIDfmincon([0.5 100 0]);
        errorPID2 = p_Zad3PID(xPID);
        
    else
    %odp skokowa dla DMC na całej przestrzeni sterowania    
        czas_sym = 600;

        u = ones(czas_sym, 1);
        y=zeros(czas_sym, 1);

        for k = 7:czas_sym
            y(k) = symulacja_obiektu3y(u(k-5),u(k-6),y(k-1),y(k-2));
        end
    %     plot(y);
    %     xlabel('k');
    %     ylabel('y(k)');
    %     hold on;
        % normalizowanie odpowiedzi skokowej
        
        st = y(2:end);
    %     regulacja DMC
        errorDMC = p_Zad3DMC([53 53 1]);
        %nastawy: [53 53 53 1] E=279,0880
        xDMC = single_DMCga();
        errorDMC2 = p_Zad3DMC(xDMC);
    end
    
else
    if DMC == 0
        %testowe argumenty
       % errorFPID = p_Zad5PIDRozm([0.5 100 0; 0.2 40 0; 0.5 100 0]);
       xPID = fuzzy_PIDfmincon2([0.5; 100; 0; 0.2; 40; 2; 0.5; 100; 2]);
       errorPID2 = p_Zad3PID(xPID);
        
    else
        st = p_odpskokFDMC(3);
        errorFDMC = p_Zad5DMCRozm(3, [53 53 1; 53 53 1; 53 53 1], st);
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