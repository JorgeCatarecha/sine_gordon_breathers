tic
km=11; % Índice máximo de los coeficientes de Fourier
wb=.8; % Frecuencia del breather
C=0; % Constante de acoplo en el Laplaciano clásico
% C0 = 1; %Constante inicial
% Ch = 0.01; %Paso en la constante fraccional
% Cff=5; % Constante final de acoplo en el Laplaciano fraccional 
m=100; % Número de partículas
%s=0.5;
n=-m/2:m/2-1;
%plist={wb,km,m,C,Cf};

z0=zeros(km+1,m); % Semilla. La solución será una matriz z(k-1,n); por tanto, la primera fila corresponde a k=0 
z0(2,m/2+1)=1;
%Si queremos empezar en una solución ya calculada
%z0 = table2array(readtable('.\Solucion_Fourier\m750s0.5km20\750C1.csv'));

options=optimoptions('fsolve','Display','iter','Jacobian','on','FunctionTolerance',1e-9,'OptimalityTolerance', 1e-9);

%Prolongamos la solución
%n_pasos = (Cff-Cf)/Ch+1;
ss = linspace(0.5,0.99,50);
for i=1:50
    s = ss(i);
    z0=zeros(km+1,m); 
    z0(2,m/2+1)=1;
    for j = linspace(0.01,2,200)
        cfrac = j;
        plist={wb,km,m,C,j};
        zk=fsolve(@(x)fbreather(x,plist,s),z0,options); %Resolvemos y usamos como semilla para la siguiente
        z0 = zk;
        % Guardar soluciones
        txt = ['.\Calculo_Bifurcacion\Soluciones\m',num2str(m),'s',num2str(s),'C',num2str(j),'.csv'];
        writematrix(zk,txt)
        % Cálculo Autovalores
        [~,D]=estabfracfourier2(100,cfrac,zk,2*pi/0.8,1,s);
        B = diag(D);
        txt = ['.\Calculo_Bifurcacion\Autovalores\m',num2str(m),'s',num2str(s),'C',num2str(j),'.csv'];
        writematrix(diag(D),txt)
        automax = max(abs(B'));
        if automax > 1.05
            break
        end
    end
        for h = linspace(0,0.01,11)
        cfrac2 = cfrac -0.01 + h;
        plist={wb,km,m,C,cfrac2};
        zk=fsolve(@(x)fbreather(x,plist,s),z0,options); %Resolvemos y usamos como semilla para la siguiente
        z0 = zk;
        % Guardar soluciones
        txt = ['.\Calculo_Bifurcacion\Soluciones\m',num2str(m),'s',num2str(s),'C',num2str(cfrac2),'.csv'];
        writematrix(zk,txt)
        % Cálculo Autovalores
        [~,D]=estabfracfourier2(100,cfrac2,zk,2*pi/0.8,1,s);
        B = diag(D);
        txt2 = ['.\Calculo_Bifurcacion\Autovalores\m',num2str(m),'s',num2str(s),'C',num2str(cfrac2),'.csv'];
        writematrix(diag(D),txt2)
        automax = max(abs(B'));
        if automax > 1.05
            cfrac2
            break
        end
        end
end


toc