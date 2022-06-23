tic
n = 101; wb = 0.8; T = 2*pi/wb;
Cf = 5;
%Número de pasos a prolongar, relacionado con donde empezamos, donde acabamos y el paso
num_pasos = 156; 
alfa = 0.5;
% Empieza desde cero si se quiere que se empiece desde una solución más
% avanzada usar soluciones en subcarpeta Soluciones
% u00 = [linspace(0,0,(n-1)/2) 1.2268 linspace(0,0,(n-1)/2)]';
% u00 = [linspace(0,0,(n-1)/2) 1.2268 linspace(0,0,(3*n-1)/2)]';
u00 = table2array(readtable('.\SolucionesTol2\N101C3.45s0.5.csv'));
u0 = u00;
Cs = linspace(3.45,Cf,num_pasos);
E = linspace(0,0,num_pasos);
mat = matfrac(n,1,alfa);
%% 
for i=1:num_pasos 
    C = Cs(i);
    u0 = contfrac2(n,wb,u0,C,mat);  
    E(i) = energfrac(u0,C,alfa);    
    % Hacer que guarde la solución en un archivo con nombre numero de nodos
    % y constante de acople.
     txt = ['.\SolucionesTol2\N',num2str(n),'C',num2str(Cs(i)),'s0.5.csv'];
     writematrix(u0,txt)
end
toc 

%% Crear animación
% animacionfrac2(u0,C,T,'FracFinal.avi',alfa)

%% Gráfica energía
% figure(2)
% plot(Cs,E,'r*')
% title('Energía')
% xlabel('alfa')
% ylabel('Energía')
