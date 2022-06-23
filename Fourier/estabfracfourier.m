function [V,D] = estabfracfourier(n,C,zk,T,bol,s)
% n: numero de nodos, C: cte acoplamiento, zk: matriz coeficientes fourier
% T: periodo = 2pi/0.8, bol = 1 (para que muestre la grafica si se quiere)
% s: fraccionalidad
tic
[coef_lim,~] = size(zk);
pos_ini=ones(1,coef_lim)*[zk(1,:);2*zk(2:end,:)];

% u es el vector de posiciones y velocidades iniciales
u = [pos_ini';zeros(n,1)];

per = [0,T];
matFlo = zeros(2*n);
m=toeplitz([0:n/2-1 n/2:-1:1]); %Problema con condiciones de contorno periódicas

% Aquí se introduce el Kernel K(n)

% Fraccional
Kernel = ker(m,s);
Kernel = Kernel - diag(diag(Kernel)); %Para quitar el término central
Kernel(isinf(Kernel))=0; %Por si el kernel diverge
mat=Kernel-sum(Kernel(1,:))*eye(n);

opcs2 = odeset('RelTol',1e-7,'AbsTol',1e-7);
for i = 1:2*n
    t = zeros(2*n,1);
    t(i) = 1;
    ini = [t;u];
    sol = ode45(@odefun,per,ini,opcs2);
    v = deval(sol,T,1:2*n);
    matFlo(i,:) = v;
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
        xi = avec(1:n)';
        dxi = avec(n+1:2*n);
        pr = sign(1i*(xi*conj(dxi)-conj(xi)*dxi));
        signat(i) = pr;
    end
    for i = 1:2*n
      if signat(i) == 1
           plot(real(d(i)),imag(d(i)),'r*')
      elseif signat(i) == -1
          plot(real(d(i)),imag(d(i)),'ro')
      else
          plot(real(d(i)),imag(d(i)),'r+')
      end
    end
    hold off
toc
end
    function [du] = odefun(~,u)
        % Ordenar primero perturbación y luego posiciones
        xi = u(1:n);
        pi1 = u(n+1:2*n);
        pos = u(2*n+1:3*n);
        vel = u(3*n+1:4*n);
        dxi = pi1;
        dpi = -cos(pos).*xi + C*mat*xi;
        dpos = vel;
        dvel = -sin(pos) + C*mat*pos;
        du = [dxi; dpi; dpos; dvel];
    end
    function [coef] = ker(m,s)
        Ls = gamma(s+0.5)*4^s./(sqrt(pi)*abs(gamma(-s)));
        Lm = gamma(abs(m)-s)./gamma(abs(m)+1+s);
        Lm(isnan(Lm)) = 0;
        coef = Ls.*Lm;
    end
end