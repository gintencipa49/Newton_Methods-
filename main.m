clear;close all;clc 
%% NOTA:LOS NODOS DEBEN SEGUIR DE MANERA ORDENADA LA SIGUIENTE
%% NUMERACION: NODO 1=SLACK(1) NODO 2..END-1=PQ(2) ULTIMO NODO=NODO PV(3)

tolerancia=1e-3;
[Ybus,Y_polar,nodos]=Ybusf(2); %el numero no afecta el resultado 
[eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos]=geneq(8); %el numero no afecta el resultado 
k=1
clc;
%% Newton completo 
disp('1 Newton Completo')
disp('2 Newton Modifcado')
disp('3 Newton Desacoplado')
disp('4 Newton Desacoplado Rapido')
n=input('Ingrese el numero del metodo: ');
switch n
    case 1
        [CI,time,iter,del_1,vartex,CI_1] = NC(eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos,tolerancia);
        clc;
        disp('Metodo de Newton_Completo ')
        disp('NOTA:LOS NODOS DEBEN SEGUIR DE MANERA ORDENADA LA SIGUIENTE NUMERACION: NODO 1=SLACK(1) NODO 2..END-1=PQ(2) ULTIMO NODO=NODO PV(3)')
        disp(['Tolerancia=',num2str(tolerancia)])
        disp(['tiempo:'])
        disp(time)
        disp(' Numero de iteraciones:')
        disp(iter)
        disp('Delta de potencia por cada iteracion (la ultima columna es la ultima iteracion) ')
        disp(del_1)
        disp('Variables en el orden de los resultados')
        disp(vartex)
        disp('Resultados por cada iteracion (la ultima columna es la ultima iteracion y los angulos estan en GRADOS):')
        [f,c]=size(CI_1);
        MIT=f-length(bui+1:bui+1+(nodos-PV-2));  
        CI_1G=CI_1;
        CI_1G(1:MIT,1:end)=CI_1(1:MIT,1:end)*180/pi;
        disp(CI_1G)
        hold on
        xlabel ('Iteraciones');
        ylabel('\Delta PQ');
        for yu=1:f
            plot(abs(del_1(yu,:)));
            yop=vartex(1:f,1:2);
            legend(eval('yop'))   ;
        end
        title('Newton Raphson Completo')
        txt = ['Tiempo: ' num2str(time) ' segundos'];
        text(2,0.5,txt)
        txt = ['Numero de iteraciones: ' num2str(iter) ''];
        text(2,0.4,txt);
        set(gcf,'color','w')
        grid on
        hold off
        disp('!IMPORTANTE¡Desplácese hacia arriba en la consola para ver los resultados')
    case 2
        [CI,time,iter,del_1,vartex,CI_1] = NM(eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos,tolerancia);
        clc;
        disp('Metodo de Newton_Completo Modificado ')
        disp('NOTA:LOS NODOS DEBEN SEGUIR DE MANERA ORDENADA LA SIGUIENTE NUMERACION: NODO 1=SLACK(1) NODO 2..END-1=PQ(2) ULTIMO NODO=NODO PV(3)')
        disp(['Tolerancia=',num2str(tolerancia)])
        disp(['tiempo:'])
        disp(time)
        disp(' Numero de iteraciones:')
        disp(iter)
        disp('Delta de potencia por cada iteracion (la ultima columna es la ultima iteracion) ')
        disp(del_1)
        disp('Variables en el orden de los resultados')
        disp(vartex)
        disp('Resultados por cada iteracion (la ultima columna es la ultima iteracion y los angulos estan en GRADOS):')
        [f,c]=size(CI_1);
        MIT=f-length(bui+1:bui+1+(nodos-PV-2));  
        CI_1G=CI_1;
        CI_1G(1:MIT,1:end)=CI_1(1:MIT,1:end)*180/pi;
        disp(CI_1G)  
        hold on
        xlabel ('Iteraciones')
        ylabel('\Delta PQ')
        for yu=1:f
            plot(abs(del_1(yu,:)));
            yop=vartex(1:f,1:2);
            legend(eval('yop'))   ;
        end
        title('Newton Raphson Modificado')
        txt = ['Tiempo: ' num2str(time) ' segundos'];
        text(2,0.5,txt);
        txt = ['Numero de iteraciones: ' num2str(iter) ''];
        text(2,0.4,txt);
        set(gcf,'color','w')
        grid on
        hold off
        disp('!IMPORTANTE¡Desplácese hacia arriba en la consola para ver los resultados')
    case 3
        [CI,time,iter,del_1,vartex,CI_1,vones] = ND(eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos,tolerancia);
        clc;
        disp('Metodo de Newton_Desacoplado')
        disp('NOTA:LOS NODOS DEBEN SEGUIR DE MANERA ORDENADA LA SIGUIENTE NUMERACION: NODO 1=SLACK(1) NODO 2..END-1=PQ(2) ULTIMO NODO=NODO PV(3)')
        disp(['Tolerancia=',num2str(tolerancia)])
        disp(['tiempo:'])
        disp(time)
        disp(' Numero de iteraciones:')
        disp(iter)
        disp('Delta de potencia por cada iteracion (la ultima columna es la ultima iteracion) ')
        disp(del_1)
        disp('Variables en el orden de los resultados')
        disp(vartex)
        disp('Resultados por cada iteracion (la ultima columna es la ultima iteracion y los angulos estan en GRADOS):')
        [f,c]=size(CI_1);
        MIT=f-length(bui+1:bui+1+(nodos-PV-2));  
        CI_1G=CI_1;
        CI_1G(1:MIT,1:end)=CI_1(1:MIT,1:end)*180/pi;
        disp(CI_1G) 
        hold on
        xlabel ('Iteraciones')
        ylabel('\Delta PQ')
        for yu=1:f
            plot(abs(del_1(yu,:)));
            yop=vartex(1:f,1:2);
            legend(eval('yop'));   
        end
        title('Newton Raphson Desacoplado')
        txt = ['Tiempo: ' num2str(time) ' segundos'];
        text(2,0.5,txt)
        txt = ['Numero de iteraciones: ' num2str(iter) ''];
        text(2.5,0.4,txt);
        set(gcf,'color','w');
        grid on
        hold off
        disp('!IMPORTANTE¡Desplácese hacia arriba en la consola para ver los resultados')
    case 4
        [CI,time,iter,del_1,vartex,CI_1] = NDR(eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos,tolerancia);
        clc;
        disp('Metodo de Newton_Desacoplado_Rapido')
        disp('NOTA:LOS NODOS DEBEN SEGUIR DE MANERA ORDENADA LA SIGUIENTE NUMERACION: NODO 1=SLACK(1) NODO 2..END-1=PQ(2) ULTIMO NODO=NODO PV(3)')
        disp(['Tolerancia=',num2str(tolerancia)])
        disp(['tiempo:'])
        disp(time)
        disp(' Numero de iteraciones:')
        disp(iter)
        disp('Delta de potencia por cada iteracion (la ultima columna es la ultima iteracion) ')
        disp(del_1)
        disp('Variables en el orden de los resultados')
        disp(vartex)
        disp('Resultados por cada iteracion (la ultima columna es la ultima iteracion y los angulos estan en GRADOS):')
        [f,c]=size(CI_1);
        MIT=f-length(bui+1:bui+1+(nodos-PV-2));  
        CI_1G=CI_1;
        CI_1G(1:MIT,1:end)=CI_1(1:MIT,1:end)*180/pi;
        disp(CI_1G) 
        hold on
        xlabel ('Iteraciones')
        ylabel('\Delta PQ')
        for yu=1:f
            plot(abs(del_1(yu,:)));
            yop=vartex(1:f,1:2);
            legend(eval('yop'));   
        end
        title('Newton Raphson Desacoplado Rapido')
        txt = ['Tiempo: ' num2str(time) ' segundos'];
        text(2.5,0.5,txt);
        txt = ['Numero de iteraciones: ' num2str(iter) ''];
        text(2.5,0.4,txt);
        set(gcf,'color','w');
        grid on
        hold off
        disp('!IMPORTANTE¡Desplácese hacia arriba en la consola para ver los resultados')
end

