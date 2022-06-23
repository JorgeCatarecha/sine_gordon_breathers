tic
n = 41; wb = 0.8; T = 2*pi/wb;
Cf = 0.001;
num_pasos = 41;
% Empieza desde cero si se quiere que se empiece desde una solución más
% avanzada usar soluciones en subcarpeta Soluciones
%u00 = [linspace(0,0,(n-1)/2) 1.2268 linspace(0,0,n-1) 0 linspace(0,0,(n-1)/2)]';
u00 = table2array(readtable('.\Soluciones\N41C0.79.csv'));
u0 = u00;
% Si se parte de una solución inicial cambiar el 0 inicial de Cs por la C
% de la solución usada
Cs = linspace(0.79,0.83,num_pasos);
E = linspace(0,0,num_pasos);

for i=1:num_pasos 
    C = Cs(i);
    u0 = cont(n,wb,u0,C);  
    E(i) = energ(u0,C);    
    % Hacer que guarde la solución, en una subcarperta Soluciones en un archivo con nombre numero de nodos
    % y constante de acople.
    txt = ['.\Soluciones\N',num2str(n),'C',num2str(Cs(i)),'.csv'];
    writematrix(u0,txt)
end
toc 

%% Crear animación
animacion(u0,C,T,'Final3.avi')

%% Gráfica energía
figure(2)
plot(Cs,E,'r*')
title('Energía')
xlabel('C')
ylabel('Energía')
