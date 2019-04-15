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
% matlab2tikz('wykresy_tikz/odpowiedzi_skok.tex', 'showInfo', false);
hold off;



ux = -1:0.1:1;
figure(2);
plot(ux,y(czas_sym, :));

% nazwa1 = sprintf('wykresy_txt/char_stat_Y.txt');
% file = fopen(nazwa1, 'w');
% B = [(1:21);y(czas_sym,:)];
% fprintf(file, '%4.3f %.3f \n',B);
% fclose(file);
% 
% nazwa1 = sprintf('wykresy_txt/char_stat_U.txt');
% file = fopen(nazwa1, 'w');
% C = [(1:21);ux];
% fprintf(file, '%4.3f %.3f \n',C);
% fclose(file);
