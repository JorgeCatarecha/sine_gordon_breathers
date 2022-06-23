function [] = animacionfrac2(u,C,T,nombre,alfa)
% Crear animación

n = length(u)/2;
mat = matfrac(n,1,alfa);
sol = ode45(@odefunfrac,[0,2*T],u);
x = linspace(1,n,n);
v = VideoWriter(nombre);
open(v);
t = linspace(0,2*T,200);
shg
for i=1:length(t)
plot(x,deval(sol,t(i),x),'r*')
ylim([-5,5])
xlabel('Nodos')
ylabel('Desplazamiento')
title(['C = ',num2str(C)])
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);

    function [du] = odefunfrac(~,u)
        pos = u(1:n);
        vel = u(n+1:2*n);
        dpos = vel;
        dvel = -sin(pos) - C*mat*pos;
        du = [dpos; dvel];
    end
end