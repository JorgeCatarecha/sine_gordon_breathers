tic
km=11; % Índice máximo de los coeficientes de Fourier
wb=.8; % Frecuencia del breather
C=0; % Constante de acoplo en el Laplaciano clásico
Cf = 0; %Constante inicial
Ch = 0.01; %Paso en la constante fraccional
Cff=1.5; % Constante final de acoplo en el Laplaciano fraccional 
m=100; % Número de partículas
s=0.5;
n=-m/2:m/2-1;
plist={wb,km,m,C,Cf};

z0=zeros(km+1,m); % Semilla. La solución será una matriz z(k-1,n); por tanto, la primera fila corresponde a k=0
 z0(2,m/2+1)=1;
%Si queremos empezar en una solución ya calculada
% z0 = table2array(readtable('.\Solucion_Fourier\m750s0.5km20\750C1.csv'));

options=optimoptions('fsolve','Display','iter','Jacobian','on','FunctionTolerance',1e-9,'OptimalityTolerance', 1e-9);

%Prolongamos la solución
n_pasos = (Cff-Cf)/Ch+1;
for i=1:n_pasos
plist{5} = (i-1)*Ch+Cf;
zk=fsolve(@(x)fbreather(x,plist,s),z0,options); %Resolvemos y usamos como semilla para la siguiente
z0 = zk;
% Guardar soluciones
 txt = ['.\Solucion_Fourier\m100s0.5km11\',num2str(m),'C',num2str((i-1)*Ch+Cf),'.csv'];
      writematrix(zk,txt)
end
t=linspace(0,2*pi/wb,51); % Valores de tiempo para la evolución temporal de la próxima línea
% El 2 que aparece multiplicando los coeficientes de Fourier, es análogo al
% a0/2 que se suele tomar (en este caso hecho al revés).
ut=(cos((0:km)'*wb*t))'*[zk(1,:);2*zk(2:end,:)]; % Evolución temporal del breather. Importante para el análisis de Floquet

%Gráfica que muestra el tiempo en eje y y el color es el desplazamiento
figure(1);imagesc(n,t,ut);axis xy
%Da la posición inicial del breather
figure(2);plot(n,ut(1,:),'.-','MarkerSize',10,'MarkerEdgeColor','r')

toc