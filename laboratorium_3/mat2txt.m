clear all;
load('../dane/pid_single.mat') 

Y = Y(1:end-1);

 figure(1)
        subplot(2,1,1);
        plot(Y);
        hold on;
        plot(yZad+Ypp);
        hold off;
        title(['Regulator PID K=',sprintf('%g',K'),' Ti=',sprintf('%g',Ti),' Td=',sprintf('%g',Td)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(U);
        
 error = sum(((yZad+Ypp) - Y).^2);

 nazwa1 = sprintf('wykresy_txt/PID_single/U__PID_K=%g_Ti=%g_Td=%g_E=%g_.txt',K,Ti,Td,error);
 nazwa2 = sprintf('wykresy_txt/PID_single/Y__PID_K=%g_Ti=%g_Td=%g_E=%g_.txt',K,Ti,Td,error);
 nazwa3 = 'wykresy_txt/PID_single/Yzad_single.txt';


 file = fopen(nazwa1, 'w');
 A = [(1:iterNum);U'];
 fprintf(file, '%4.3f %.3f \n',A);
 fclose(file);

 file = fopen(nazwa2, 'w');
 B = [(1:iterNum);Y'];
 fprintf(file, '%4.3f %.3f \n',B);
 fclose(file);

 file = fopen(nazwa3, 'w');
 C = [(1:iterNum);(yZad+Ypp)'];
 fprintf(file, '%4.3f %.3f \n',C);
 fclose(file);
 
  %----------------------------------------------------------------
  clear all;
 load('../dane/dmc_single.mat') 
 
 Y = Y(1:end-1);

 figure(1)
        subplot(2,1,1);
        plot(Y);
        hold on;
        plot(yZad+Ypp);
        hold off;
        title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(U);
        
 error = sum(((yZad+Ypp) - Y).^2);

 nazwa1 = sprintf('wykresy_txt/DMC_single/U__DMC_D=%g_N=%g_Nu=%g_L=%g_E=%g_.txt',D,N,Nu,lambda,error);
 nazwa2 = sprintf('wykresy_txt/DMC_single/Y__DMC_D=%g_N=%g_Nu=%g_L=%g_E=%g_.txt',D,N,Nu,lambda,error);
 nazwa3 = 'wykresy_txt/DMC_single/Yzad_single.txt';


 file = fopen(nazwa1, 'w');
 A = [(1:iterNum);U'];
 fprintf(file, '%4.3f %.3f \n',A);
 fclose(file);

 file = fopen(nazwa2, 'w');
 B = [(1:iterNum);Y'];
 fprintf(file, '%4.3f %.3f \n',B);
 fclose(file);

 file = fopen(nazwa3, 'w');
 C = [(1:iterNum);(yZad+Ypp)'];
 fprintf(file, '%4.3f %.3f \n',C);
 fclose(file);
 
 %----------------------------------------------------------------
 clear all;
 load('../dane/fuzzy_pid_k=8,8,7_Ti=60,50,40.mat') 

 Y = Y(1:end-1);
 figure(1)
        subplot(2,1,1);
        plot(Y);
        hold on;
        plot(yZad+Ypp);
        hold off;
        title(['Regulator PID K=',sprintf('%g',K'),' Ti=',sprintf('%g',Ti),' Td=',sprintf('%g',Td)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(Ukonc);
        
 error = sum(((yZad+Ypp) - Y).^2);

 nazwa1 = sprintf('wykresy_txt/PID_fuzzy/U__PID_K=%g,%g,%g_Ti=%g,%g,%g_Td=%g,%g,%g_E=%g_.txt',K,Ti,Td,error);
 nazwa2 = sprintf('wykresy_txt/PID_fuzzy/Y__PID_K=%g,%g,%g_Ti=%g,%g,%g_Td=%g,%g,%g_E=%g_.txt',K,Ti,Td,error);
 nazwa3 = 'wykresy_txt/PID_fuzzy/Yzad_fuzzy.txt';


 file = fopen(nazwa1, 'w');
 A = [(1:iterNum);Ukonc'];
 fprintf(file, '%4.3f %.3f \n',A);
 fclose(file);

 file = fopen(nazwa2, 'w');
 B = [(1:iterNum);Y'];
 fprintf(file, '%4.3f %.3f \n',B);
 fclose(file);

 file = fopen(nazwa3, 'w');
 C = [(1:iterNum);(yZad+Ypp)'];
 fprintf(file, '%4.3f %.3f \n',C);
 fclose(file);
 
  %----------------------------------------------------------------
  clear all;
 load('../dane/fuzzy_pid_k=9,8,7_Ti=50,70,60.mat') 
 
 Y = Y(1:end-1);
 figure(2)
        subplot(2,1,1);
        plot(Y);
        hold on;
        plot(yZad+Ypp);
        hold off;
        title(['Regulator PID K=',sprintf('%g',K'),' Ti=',sprintf('%g',Ti),' Td=',sprintf('%g',Td)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(Ukonc);
        
 error = sum(((yZad+Ypp) - Y).^2); 

 nazwa1 = sprintf('wykresy_txt/PID_fuzzy/U__PID_K=%g,%g,%g_Ti=%g,%g,%g_Td=%g,%g,%g_E=%g_.txt',K,Ti,Td,error);
 nazwa2 = sprintf('wykresy_txt/PID_fuzzy/Y__PID_K=%g,%g,%g_Ti=%g,%g,%g_Td=%g,%g,%g_E=%g_.txt',K,Ti,Td,error);
 nazwa3 = 'wykresy_txt/PID_fuzzy/Yzad_fuzzy.txt';


 file = fopen(nazwa1, 'w');
 A = [(1:iterNum);Ukonc'];
 fprintf(file, '%4.3f %.3f \n',A);
 fclose(file);

 file = fopen(nazwa2, 'w');
 B = [(1:iterNum);Y'];
 fprintf(file, '%4.3f %.3f \n',B);
 fclose(file);

 file = fopen(nazwa3, 'w');
 C = [(1:iterNum);(yZad+Ypp)'];
 fprintf(file, '%4.3f %.3f \n',C);
 fclose(file);
 
 %----------------------------------------------------------------
 clear all;
 load('../dane/fuzzy_dmc_lambda=1,1,1.mat') 

 Y = Y(1:end-1);
 figure(3)
        subplot(2,1,1);
        plot(Y);
        hold on;
        plot(yZad+Ypp);
        hold off;
        title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(Ukonc);
        
 error = sum(((yZad+Ypp) - Y).^2);

 nazwa1 = sprintf('wykresy_txt/DMC_fuzzy/U__DMC_D=%g_N=%g,%g,%g_Nu=%g,%g,%g_L=%g,%g,%g_E=%g_.txt',D',N',Nu',lambda',error);
 nazwa2 = sprintf('wykresy_txt/DMC_fuzzy/Y__DMC_D=%g_N=%g,%g,%g_Nu=%g,%g,%g_L=%g,%g,%g_E=%g_.txt',D',N',Nu',lambda',error);
 nazwa3 = 'wykresy_txt/DMC_fuzzy/Yzad_fuzzy.txt';


 file = fopen(nazwa1, 'w');
 A = [(1:iterNum);Ukonc];
 fprintf(file, '%4.3f %.3f \n',A);
 fclose(file);

 file = fopen(nazwa2, 'w');
 B = [(1:iterNum);Y'];
 fprintf(file, '%4.3f %.3f \n',B);
 fclose(file);

 file = fopen(nazwa3, 'w');
 C = [(1:iterNum);(yZad+Ypp)'];
 fprintf(file, '%4.3f %.3f \n',C);
 fclose(file);
 
 %----------------------------------------------------------------
 clear all;
 load('../dane/fuzzy_dmc_lambda=1,10,8.mat') 

 Y = Y(1:end-1);
 figure(4)
        subplot(2,1,1);
        plot(Y);
        hold on;
        plot(yZad+Ypp);
        hold off;
        title(['Regulator DMC D=',sprintf('%g',D'),' N=',sprintf('%g',N),' Nu=',sprintf('%g',Nu),' lambda=',sprintf('%g',lambda)]);
        legend('y','yzad')
        subplot(2,1,2);
        stairs(Ukonc);
        
 error = sum(((yZad+Ypp) - Y).^2);

 nazwa1 = sprintf('wykresy_txt/DMC_fuzzy/U__DMC_D=%g_N=%g,%g,%g_Nu=%g,%g,%g_L=%g,%g,%g_E=%g_.txt',D',N',Nu',lambda',error);
 nazwa2 = sprintf('wykresy_txt/DMC_fuzzy/Y__DMC_D=%g_N=%g,%g,%g_Nu=%g,%g,%g_L=%g,%g,%g_E=%g_.txt',D',N',Nu',lambda',error);
 nazwa3 = 'wykresy_txt/DMC_fuzzy/Yzad_fuzzy.txt';

 file = fopen(nazwa1, 'w');
 A = [(1:iterNum);Ukonc];
 fprintf(file, '%4.3f %.3f \n',A);
 fclose(file);

 file = fopen(nazwa2, 'w');
 B = [(1:iterNum);Y'];
 fprintf(file, '%4.3f %.3f \n',B);
 fclose(file);

 file = fopen(nazwa3, 'w');
 C = [(1:iterNum);(yZad+Ypp)'];
 fprintf(file, '%4.3f %.3f \n',C);
 fclose(file);
 