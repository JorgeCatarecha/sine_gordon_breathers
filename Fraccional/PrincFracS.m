tic
% Continúa la solución clásica al pasar al caso fraccional
n = 101; wb = 0.8; T = 2*pi/wb;
C = 1.5;
num_pasos = 23;
alfaf = 0.5;
% Cálculo solución inicial a partir solución clásica
f = @(x) 4*atan(sqrt(1-wb^2)./(wb*cosh(sqrt(1-wb^2)*x)));
% u00 = zeros(2*n,1);
% dist = sqrt(1/C);
% for i=1:n
%    u00(i) = f((i-(n+1)/2)*dist);
% end
% Solución inicial ya calculada
u00 = table2array(readtable('.\SolucionesConts\N101C1.5s0.72.csv'));
%%
u0= u00;
alfas = linspace(0.72,alfaf,num_pasos);
E = linspace(0,0,num_pasos);

for i=1:num_pasos 
    alfa = alfas(i);
    u0 = conts(n,wb,u0,C,alfa);  
    E(i) = energfrac(u0,C,alfa);    
    % Hacer que guarde la solución en un archivo con nombre numero de nodos
    % y constante de acople.
     txt = ['.\SolucionesConts\N',num2str(n),'C',num2str(C),'s',num2str(alfa),'.csv'];
     writematrix(u0,txt)
end
toc 

%% Crear animación
animacionfrac2(u0,C,T,'FracFinalS.avi',alfa)

%% Gráfica energía
figure(2)
plot(alfas,E,'r*')
title('Energía')
xlabel('Alfa')
ylabel('Energía')
