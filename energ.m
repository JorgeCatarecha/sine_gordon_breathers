function [E] = energ(u,C)
% Calcula la energía de la solución u, con una constante de acoplo C
n = length(u)/2;
pos = u(1:n); 
vel = u(n+1:2*n);
Ecin = sum(vel.^2)/2;
U1 = sum((pos-circshift(pos,1)).^2)/2; %Energía potencial elástica, caso periódico
U2 = sum(1-cos(pos)); %Energía potencial tipo coseno
E = Ecin + C*U1 + U2;
end