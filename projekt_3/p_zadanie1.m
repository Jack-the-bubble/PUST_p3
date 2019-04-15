% Zadanko pierwsze

czas_sym = 30;


u= zeros(czas_sym, 1);
y=zeros(4,czas_sym);
for i = 0:3
    y(i+1,:)=i;
    for k = 7:czas_sym
        y(i+1,k) = symulacja_obiektu3y(u(k-5),u(k-6),y(i+1,k-1),y(i+1,k-2));
    end

end

figure(1)
plot(y');
xlabel('k');
ylabel('y(k)');
hold on;
matlab2tikz('wykresy_tikz/punkt_pracy.tex', 'showInfo', false);

% file = fopen(nazwa3, 'w');
% C = [(1:czas_sym);yzad'];
% fprintf(file, '%4.3f %.3f \n',C);
% fclose(file);
