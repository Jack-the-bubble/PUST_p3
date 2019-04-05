clear all;
DMC = 0;

if DMC == 0
    p_Zad3PID([0.5 100 0]);


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
    p_Zad3DMC([53 53 1], st);
    
end