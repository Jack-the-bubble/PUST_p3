%PUST Lab3
%Charakterystyka statyczna

u = [20 30 40 50 60 70 80];
y = [28.31 33 37.93 42.56 45.25 47.56 50.06];

plot(u,y,'o-');
matlab2tikz('wykresy_tikz/char_stat.tex', 'showInfo', false);
 
% nazwa1 = 'wykresy_txt/char_stat.txt';
% 
%  file = fopen(nazwa1, 'w');
%  A = [u;y];
%  fprintf(file, '%4.3f %.3f \n',A);
%  fclose(file);