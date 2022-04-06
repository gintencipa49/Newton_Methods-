function [eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos] = geneq(x)
%Matlab R2021b academic use 
%GENEQ Summary of this function goes here
%   Funcion que permite obtener las ecuaciones de cualquier sistema de
%   potencia (los datos deben estar ordenados donde el primero sea el nodo
%   slack y los ultimos los nodos PQ
[Ybus]=Ybusf(2); 
x=readtable('datos_nodos.txt');
X=table2array(x);
x_1=readtable('datos_linea.txt');               %PEDIR AL USUARIO 
X1=table2array(x_1);
a=size(X);
nodo=max(X(:,1:2));
nodos=max(nodo);
tipoo=find(X(:,2)==3);
PV=length(tipoo);
%% varibales simbolicas
vi=sym('v%d',[1,a(1,1)]);
oi=sym('o%d',[1,a(1,1)]);
pit=sym('p%d',[1,a(1,1)]);
qi=sym('q%d',[1,a(1,1)]);
%% lectura de datos 
% potencias 
p=zeros(a(1,1),1);
q=zeros(a(1,1),1);
v=ones(a(1,1),1);
o=zeros(a(1,1),1);
vf=zeros(a(1,1),1);
of=zeros(a(1,1),1);
for ju=1:a(1,1)
    ki=X(ju,2);
    if ki==1
        vi(1,ju)=X(ju,3);
        oi(1,ju)=X(ju,4); 
    end 
    if ki==2
        %PQ
        pit(1,ju)=-X(ju,6);
        qi(1,ju)=-X(ju,7);
    end
    if ki==3
        %PV
        pit(1,ju)=X(ju,5);
        vi(1,ju)=X(ju,3);
    end 
end
% p coienza en 2 y termina en el ultimo nodo 
% q comienza en 2 y y termina en el ultimo nodo-1
%% potencias esperadas 
PQ=pit(1,2:end);
bui=length(PQ);
PQ(1,bui+1:bui+1+(nodos-PV-2))=qi(1,2:end-PV);
%% condiciones iniciales 
% v coienza en 2 y termina en el ultimo nodo 
% o comienza en 2 y y termina en el ultimo nodo-1
CI=o(2:end,1);
CI(bui+1:bui+1+(nodos-PV-2),1)=v(2:end-1,1);
iniciales=CI;
%% DEFINIR VARIABLES SYMBOLICAS 
%text='syms '
for hy=2:nodos      %Nodo 1 siempre debe ser slack
    ecuacion=0;
    busc=find(X1(:,1:2)==hy);
    con=length(busc);
    [fila,columna]=find(X1==2);
    cor=[fila columna];
    conex=X1(fila, 1:2);
    for lh=1:con
            [fila,columna]=find(X1==hy);
            cor=[fila columna];
            %conex(oi-1,2)=X1(cor(oi,1),cor(oi,2));
             if cor(lh,2)==1
                conex(lh,2)=X1(cor(lh,1),cor(lh,2)+1); 
             end
             if cor(lh,2)==2
                conex(lh,2)=X1(cor(lh,1),cor(lh,2)-1);
             end
             conex(lh,1)=hy;
    end
    for lm=1:con 
        eq=ecuacion+(abs(Ybus(conex(lm,1),conex(lm,2))*vi(1,conex(lm,1))*vi(1,conex(lm,2)))*...
            cos( oi(1,conex(lm,1)) - oi(1,conex(lm,2))-angle(Ybus(conex(lm,1),conex(lm,2)))));
        ecuacion=eq;
    end
    eq=eq+(abs((vi(1,conex(lm,1)))^2*Ybus(conex(lm,1),conex(lm,1))))*cos(-angle(Ybus(conex(lm,1),conex(lm,1))));
    eqt(hy-1,1)=subs(eq);
end
memoria=hy;
for hy=2:nodos-PV      %Nodo 1 siempre debe ser slack
    ecuacion=0;
    busc=find(X1(:,1:2)==hy);
    con=length(busc);
    [fila,columna]=find(X1==2);
    cor=[fila columna];
    conex=X1(fila, 1:2);
    for lh=1:con
            [fila,columna]=find(X1==hy);
            cor=[fila columna];
            %conex(oi-1,2)=X1(cor(oi,1),cor(oi,2));
             if cor(lh,2)==1
                conex(lh,2)=X1(cor(lh,1),cor(lh,2)+1); 
             end
             if cor(lh,2)==2
                conex(lh,2)=X1(cor(lh,1),cor(lh,2)-1);
             end
             conex(lh,1)=hy;
    end
    for lm=1:con 
        eq=ecuacion+(abs(Ybus(conex(lm,1),conex(lm,2))*vi(1,conex(lm,1))*vi(1,conex(lm,2)))*...
            sin( oi(1,conex(lm,1)) - oi(1,conex(lm,2))-angle(Ybus(conex(lm,1),conex(lm,2)))));
        ecuacion=eq;
    end
    eq=eq+(abs((vi(1,conex(lm,1)))^2*Ybus(conex(lm,1),conex(lm,1))))*sin(-angle(Ybus(conex(lm,1),conex(lm,1))));
    eqt(memoria,1)=subs(eq);
    memoria=memoria+1;
end
