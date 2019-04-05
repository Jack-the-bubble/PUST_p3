clear all;
czas_sym = 200;


u= ones(czas_sym, 1);
y=zeros(czas_sym, 21);
figure(1);
for i = 1:21
    u=ones(czas_sym, 1)*(-1.1+i*0.1);
    
    for k = 7:czas_sym
        y(k, i) = symulacja_obiektu3y(u(k-5),u(k-6),y(k-1, i),y(k-2, i));
    end
    plot(y);
    xlabel('k');
    ylabel('y(k)');
    hold on;
end
hold off;

ux = -1:0.1:1;

figure(2);
plot(ux,y(czas_sym, :));

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
