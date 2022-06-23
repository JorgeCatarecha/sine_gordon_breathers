function [] = animacion(u,C,T,nombre)
% Crear animación
n = length(u)/2;
sol = ode45(@odefun,[0,2*T],u);
x = linspace(1,n,n);
v = VideoWriter(nombre);
open(v);
t = linspace(0,2*T,200);
shg
for i=1:length(t)
plot(x,deval(sol,t(i),x),'r*')
ylim([-3,3])
xlabel('Nodos')
ylabel('Desplazamiento')
title(['C = ',num2str(C)])
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);

    function [du] = odefun(~,u)
        pos = u(1:n);
        vel = u(n+1:2*n);
        dpos = vel;
        dvel = -sin(pos) - C*(2*pos-circshift(pos,1)-circshift(pos,-1));
        du = [dpos; dvel];
    end
end