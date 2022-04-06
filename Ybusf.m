function [Ybus,Y_polar,l] = Ybusf(x)
%YBUS Summary of this function goes here
%   Funcion que permite obtener la matriz de un sistema de potencia, los
%   datos de entrada deben ser especificados en el archivo datos de linea;
x=readtable('datos_linea.txt');  %datos txt 
X=table2array(x);
l=max(X(:,1:2)); %Nodos
l=max(l);
Ybus=zeros(l);
k=0;
uh=[];
bk=0;
Ybus(eye(size(Ybus))==1)=1;
r=0;
%Diagonal
while k< l
    k=k+1;
    r=r+1;
    [fila,columna]=find(X==k);
    jk=length(fila);
    for i=1:jk
        uh(k)=(((X(fila(i),3))+(X(fila(i),4)))^-1)+(X(fila(i),5));
        bk=bk+uh(k);
        Ybus(r,r)=bk;
    end
for g=1:l
end
   uh=uh*0;
   bk=0;

end
% No diagonal 
for lo=1:length(X)
    kip=X(lo,1);
    jup=X(lo,2);
    Ybus(kip,jup)=-1*(((X(lo,3))+(X(lo,4)))^-1);
    Ybus(jup,kip)=Ybus(kip,jup);
end

%% Carteciano a polar 
                %Ib=196.8289*exp(i*-83.130*pi/180);
Y_mag=abs(Ybus);
Y_ang=angle(Ybus);
Y_polar=zeros(length(Ybus),2*length(Ybus));
for ty=1:length(Y_polar)
    for hi=1:length(Y_polar)/2
        if mod(ty,2)~=0
            Y_polar(hi,ty)=Y_mag(hi,(ty+1)/2) ;
        end  
        if mod(ty,2)==0
            Y_polar(hi,ty)=Y_ang(hi,ty/2)*180/pi ;
        end  
    end 
end
end

