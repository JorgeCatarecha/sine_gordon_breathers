function [E] = energfrac(u,C,alfa)
% Calcula la energía de la solución u, con una constante de acoplo C
% cc Dirichlet
n = length(u)/2;
pos = u(1:n); 
vel = u(n+1:2*n);
Ecin = sum(vel.^2)/2;
U1 = 0;
for i=1:n
U1 = U1 + energiapot(i); %Energía potencial elástica
end
U2 = sum(1-cos(pos)); %Energía potencial tipo coseno
E = Ecin + C*U1/2 + U2;
    function [energia] = energiapot(i)
        energia = 0;
        for j=i+1:n
            energia = energia + ker(abs(j-i),alfa)*(pos(i)-pos(j)).^2;
        end
    end
end