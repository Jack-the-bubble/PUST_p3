function st = p_odpskokFDMC(lreg)
%odp skokowe dla DMC na poszczegolnych kawalkach sterowania
    
    czas_sym = 600;
    chwila_skoku = 400;
    st = zeros(lreg,czas_sym-chwila_skoku);
    uzad = ones(lreg,2);
    
    
    if lreg == 2
        uzad = [];
    elseif lreg == 3
        uzad = [-0.5 -0.495; 0.2 0.205; 0.8 0.805];
    else
        disp("NIE");
        return;
    end
    
    for i=1:lreg
        u = ones(czas_sym, 1)*uzad(i,1); %stabilizacja wartosci poczatkowej
        u(chwila_skoku:end) = uzad(i,2); %skok o 0.005
        y = zeros(czas_sym, 1);

        for k = 7:czas_sym
            y(k) = symulacja_obiektu3y(u(k-5),u(k-6),y(k-1),y(k-2));
        end

        % normalizowanie odpowiedzi skokowej
        st(i,:) = (y(chwila_skoku+1:end) - y(chwila_skoku-1)) ./ (uzad(i,2) - uzad(i,1));
    end

%     figure(5)
%     for es=1:lreg
%         plot(st(es,:));
%         hold on;
%     end
%     hold off;
    
end