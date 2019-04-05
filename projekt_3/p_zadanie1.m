% Zadanko pierwsze

czas_sym = 30;


u= zeros(czas_sym, 1);
y=ones(czas_sym, 1)*10;
for i = 0:3
    y=ones(czas_sym, 1)*i;
    for k = 7:czas_sym
        y(k) = symulacja_obiektu3y(u(k-5),u(k-6),y(k-1),y(k-2));
    end
    plot(y);
    xlabel('k');
    ylabel('y(k)');
    hold on;
end
hold off;

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
