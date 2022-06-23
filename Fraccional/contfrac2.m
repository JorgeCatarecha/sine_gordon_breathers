function [uf] = contfrac2(n,wb,u,C,mat)
% Argumentos: número de nodos, frecuencia breather, y punto inicial para
% buscar el breather con constante de acople C.

T = 2*pi/wb;
per = [0,T];
opciones = optimset('Display','iter','Algorithm','levenberg-marquardt'); 
uf = fsolve(@minf,u,opciones); 
    function [dif] = minf(x0)
        opcs2 = odeset('RelTol',1e-8,'AbsTol',1e-8);
        sol = ode45(@odefunfrac,per,x0,opcs2);
        dif = norm(deval(sol,T)-x0);
   end

    function [du] = odefunfrac(~,u)
        pos = u(1:n);
        vel = u(n+1:2*n);
        dpos = vel;
        dvel = -sin(pos) - C*mat*pos;
        du = [dpos; dvel];
    end
end