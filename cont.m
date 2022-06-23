function [uf] = cont(n,wb,u0,C)
% Argumentos: número de nodos, frecuencia breather, y punto inicial para
% buscar el breather con constante de acople C.
T = 2*pi/wb;
per = [0,T];
opciones = optimset('Display','iter','TolFun',1e-8); 
uf = fsolve(@minf,u0,opciones); 
    function [dif] = minf(x0)
        opcs2 = odeset('RelTol',1e-8,'AbsTol',1e-8);
        sol = ode45(@odefun,per,x0,opcs2);
        dif = x0 - deval(sol,T);
    end

    function [du] = odefun(~,u)
        n = length(u)/2;
        pos = u(1:n);
        vel = u(n+1:2*n);
        dpos = vel;
        dvel = -sin(pos) - C*(2*pos-circshift(pos,1)-circshift(pos,-1));
        du = [dpos; dvel];
    end
end