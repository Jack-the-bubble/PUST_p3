%projekt 2xreg
u = (-1:0.01:1);
value = trapmf(u,[-1 -1 0 0.6]);
value1 = trapmf(u,[0 0.6 1 1]);
figure(1);
hold on;
plot(u,value);
plot(u,value1);
hold off;

nazwa1 = 'projekt_2xreg_1.txt';
nazwa2 = 'projekt_2xreg_2.txt';

file = fopen(nazwa1, 'w');
A = [u;value];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen(nazwa2, 'w');
B = [u;value1];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

%projekt 3xreg
value = trapmf(u,[-1 -1 -0.2 0]);
value1 = trapmf(u,[-0.2 0 0.4 0.6]);
value2 = trapmf(u,[0.4 0.6 1 1]);
figure(2);
hold on;
plot(u,value);
plot(u,value1);
plot(u,value2);
hold off;

nazwa1 = 'projekt_3xreg_1.txt';
nazwa2 = 'projekt_3xreg_2.txt';
nazwa3 = 'projekt_3xreg_3.txt';

file = fopen(nazwa1, 'w');
A = [u;value];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen(nazwa2, 'w');
B = [u;value1];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

file = fopen(nazwa3, 'w');
C = [u;value2];
fprintf(file, '%4.3f %.3f \n',C);
fclose(file);
 
%laby 3xreg
u = (0:0.1:100)
value = trapmf(u,[0 0 40 45]);
value1 = trapmf(u,[40 45 55 60]);
value2 = trapmf(u,[55 60 100 100]);
figure(3);
hold on;
plot(u,value);
plot(u,value1);
plot(u,value2);
hold off;

nazwa1 = 'lab_3xreg_1.txt';
nazwa2 = 'lab_3xreg_2.txt';
nazwa3 = 'lab_3xreg_3.txt';

file = fopen(nazwa1, 'w');
A = [u;value];
fprintf(file, '%4.3f %.3f \n',A);
fclose(file);

file = fopen(nazwa2, 'w');
B = [u;value1];
fprintf(file, '%4.3f %.3f \n',B);
fclose(file);

file = fopen(nazwa3, 'w');
C = [u;value2];
fprintf(file, '%4.3f %.3f \n',C);
fclose(file);