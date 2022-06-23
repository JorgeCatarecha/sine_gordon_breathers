function [V,D] = estab(n,C,u,T,bol)
% La función u inicial para resolver, debe ser el breather para el C acorde
per = [0,T];
matFlo = zeros(2*n);
opcs2 = odeset('RelTol',1e-10,'AbsTol',1e-10);
for i = 1:2*n
    t = zeros(2*n,1);
    t(i) = 1;
    ini = [t;u];
    sol = ode45(@odefun,per,ini,opcs2);
    v = deval(sol,T,1:2*n);
    matFlo(:,i) = v;
end

[V,D] = eig(matFlo);
if bol == 1
   % Los representa en el circulo unidad 
    d = diag(D);
    z = linspace(0,2*pi);
    plot(cos(z),sin(z),'k')
    axis equal
    hold on
    xlabel('Parte Real','FontSize',15)
    ylabel('Parte Imaginaria','FontSize',15)
    xtickformat('%.1f')
    ytickformat('%.1f')
    title('Autovalores Operador de Floquet','FontSize',15)
    grid on
    grid minor
    % Cálculo de la signatura de Krein
    signat = zeros(1,2*n);
    for i = 1:2*n
        avec = V(:,i);
        xi = avec(1:n);
        dxi = avec(n+1:2*n);
       simpl = 1i*(dxi'*xi-xi'*dxi);
        %simpl = 1i*(xi'*conj(dxi)-conj(xi')*dxi);
        pr = sign(simpl);
        signat(i) = pr;
        if abs(simpl) < 1e-3
            signat(i) = 0;
        end
    end
    for i = 1:2*n
      if signat(i) == 1
           plot(real(d(i)),imag(d(i)),'b*')
      elseif signat(i) == -1
          plot(real(d(i)),imag(d(i)),'ro')
      else
          plot(real(d(i)),imag(d(i)),'k+')
      end
    end
    hold off
end
    function [du] = odefun(~,u)
        % Ordenar primero perturbación y luego posiciones
        xi = u(1:n);
        pi1 = u(n+1:2*n);
        pos = u(2*n+1:3*n);
        vel = u(3*n+1:4*n);
        dxi = pi1;
        dpi = -cos(pos).*xi - C*(2*xi-circshift(xi,1)-circshift(xi,-1));
        dpos = vel;
        dvel = -sin(pos) - C*(2*pos-circshift(pos,1)-circshift(pos,-1));
        du = [dxi; dpi; dpos; dvel];
    end
end