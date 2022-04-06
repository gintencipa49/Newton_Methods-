function [CI,time,iter,del_1,vartex,CI_1,vones] = ND(eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos,tolerancia)
variables=oi(1,2:end);
variables(1,bui+1:bui+1+(nodos-PV-2))=vi(1,2:end-PV);
for ty=1:length(variables)
    variablest(ty,1)=variables(1,ty);
end
%% convertir a texto 
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
Jacobiano=jacobian(eqt,variables);  
tic
while error>=tolerancia
        CItex=num2str(CI,16);
        for pir=1:length(vartex)
            eval([vartex(pir,1:end) iqualtex(pir,1:end) CItex(pir,1:end)]);
        end     
        deltaPQ=double(PQ'-eval(eqt));
        del_1(:,iter+1)=double(deltaPQ);   %%opcional
        jacon=eval(Jacobiano);
        Jacov=eval(Jacobiano);
        N=nodos;
        Jacov(N:end,1:N-1)=0;   %J
        Jacov(1:N-1,N:end)=0;   %N
        Jaco=Jacov;
        Jacf=Jaco(:,nodos:end);
        [f,c]=size(Jacf);
        for pou=1:c
            Jaci(:,pou)=Jacf(:,pou).*CI(nodos-1+pou,1);
        end
        Jaco(:,nodos:end)=Jaci;
        Jac=Jaco;
        vones=ones(length(CI),1);
        vones(nodos:end,1)=CI(nodos:end,1);
        CI=CI+inv(Jac)*deltaPQ.*vones;
        error=max(abs(deltaPQ));
        CI_1(:,iter+1)=double(CI);   %%opcional
        iter=iter+1;
end
time=toc;
end

