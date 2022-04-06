function [CI,time,iter,del_1,vartex,CI_1] = NDR(eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos,tolerancia)
%NDR Summary of this function goes here
%   Detailed explanation goes here
[Ybus,Y_polar,nodos]=Ybusf(2); %el numero no afecta el resultado 
variables=oi(1,2:end);
variables(1,bui+1:bui+1+(nodos-PV-2))=vi(1,2:end-PV);
for ty=1:length(variables)
    variablest(ty,1)=variables(1,ty);
end
%% convertir a texto 
%
%iqualtex=['=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'=';'='];
for jud=1:nodos*3;
    iqualtex(jud,1)='=';
end
texxx=cellstr(string(variablest));
vartex=cell2mat(texxx);

%% Newton completo 
%% Newton completo 
error=100;
iter=0;
N=nodos;
NPQ=N-PV-1
H=-imag(Ybus(2:N,2:N));
L=-imag(Ybus(2:end-PV,2:end-PV));
N1=zeros(N-1,NPQ);
J=zeros(NPQ,N-1);
Jacobiano=[H N1; J L];
  
%%
tic
%for ddddd=1:1
while error>=tolerancia
        %%generar vectos vvvvvv
        CI_PV=double(vi(1,length(vi)-PV+1:end))';
        CI_f=[CI;CI_PV];
        kiu=CI_f(N:end,1);
        %for piu=1:length(CI)
        piu=1;
        cont=0;
        while piu<=length(CI_f)-PV
            cont=cont+1;
            vis(piu,1)=kiu(cont,1);
            if cont==NPQ+1
               cont=0;
            end
            piu=piu+1;
        end  
        %%%%%%%%%%%%%%%%%%%%%%%%
        CItex=num2str(CI,16);
        for pir=1:length(vartex)
            eval([vartex(pir,1:end) iqualtex(pir,1:end) CItex(pir,1:end)]);
        end     
        deltaPQ=double(PQ'-eval(eqt));

        DeltaPQ_V=deltaPQ./vis;
        del_1(:,iter+1)=double(DeltaPQ_V);   %%opcional
        jacon=Jacobiano;
        Jacov=Jacobiano;
        Jac=Jacov;
        CI=CI+inv(Jac)*DeltaPQ_V;
        error=max(abs(DeltaPQ_V));
        CI_1(:,iter+1)=double(CI);   %%opcional
        iter=iter+1;

end
time=toc;

end

