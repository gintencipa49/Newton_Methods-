function [CI,time,iter,del_1,vartex,CI_1] = NC(eqt,vi,oi,pit,qi,PQ,iniciales,CI,PV,bui,nodos,tolerancia)
%Definir variables cos
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
tic
while error>=tolerancia
        CItex=num2str(CI,16);
        for pir=1:length(vartex)
            eval([vartex(pir,1:end) iqualtex(pir,1:end) CItex(pir,1:end)]);
        end
        Jacobiano=jacobian(eqt,variables);
        deltaPQ=double(PQ'-eval(eqt));
        del_1(:,iter+1)=double(deltaPQ);   %%opcional
        Jac=eval(Jacobiano);
        CI=CI+inv(Jac)*deltaPQ;
        error=max(abs(deltaPQ));
        CI_1(:,iter+1)=double(CI);   %%opcional
        iter=iter+1;
end
time=toc
end

