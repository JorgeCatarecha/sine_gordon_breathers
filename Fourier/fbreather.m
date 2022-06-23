function [f,Jac]=fbreather(z,plist,s)

[wb,km,m,C,Cf]=plist{1,:};

z=reshape(z,km+1,m);

Lap=spdiags([ones(m,2) -2*ones(m,1) ones(m,2)],[-m+1 -1 0 1 m-1],m,m);

n=toeplitz([0:m/2-1 m/2:-1:1]); %Problema con condiciones de contorno periódicas

% Aquí se introduce el Kernel K(n)
% Kernel=1./abs(n).^3; % Kernel prueba
% Fraccional
Kernel = ker(n,s);
Kernel = Kernel - diag(diag(Kernel)); %Para quitar el término central
% Continuación
Kernel(isinf(Kernel))=0;
FLap=Kernel-sum(Kernel(1,:))*eye(m);

Vprimek=zeros(2*km+1,m);
for i=1:m
    z0=zeros(2*km+1,1);
    z0(1:km+1)=z(:,i);
    u=idftc(z0);
    vprime=sin(u); % Aquí viene la primera derivada del potencial
    Vprimek(:,i)=dftc(vprime);
end

f=-wb^2*(diag((0:km).^2)*z) - C*z*Lap - Cf*z*FLap + Vprimek(1:km+1,:);

f=f(:);

Jac=jac_oscm(z,plist,Lap,FLap);

end

function Jac=jac_oscm(z,plist,Lap,FLap)

[wb,km,m,C,Cf]=plist{1,:};

kk=km+1;

z=reshape(z,kk,m);
Jac=spalloc(kk*m,kk*m,(kk+m-1)*kk*m);

for n=1:m
    Jac((n-1)*kk+1:n*kk,(n-1)*kk+1:n*kk)=jac_osc1(z(:,n),wb);
end

Jac=Jac-kron(C*Lap+Cf*FLap,speye(kk));

end

function J=jac_osc1(z,wb)

km=length(z)-1;

J=zeros(km+1);
z0=zeros(2*(km+1),1);
z0(1:km+1)=z;
x=idftc(z0);
vsec=cos(x); % Aquí viene la segunda derivada del potencial
Vseck=dftc(vsec);

for k=0:km
  J(k+1,1)=Vseck(k+1);
  for l=1:km
    J(k+1,l+1)=Vseck(abs(k-l)+1)+Vseck(k+l+1);
  end
end
J=J-diag((0:km).^2)*wb^2;

end

function cm=dftc(xk)
xk=xk(:);
cm=xk;
km=length(xk)-1;
n=2*km+1;
mm=1:km;
kk=mm';
bkm=2*cos(2*pi/n*kk*mm);
xk2=xk(2:km+1);
cm(1)=(xk(1)+2*sum(xk2))/n;
cm2=1/n*(xk(1)+bkm*xk2);
cm=[cm(1);cm2];
end
function [coef] = ker(m,s)
Ls = gamma(s+0.5)*4^s./(sqrt(pi)*abs(gamma(-s)));
Lm = gamma(abs(m)-s)./gamma(abs(m)+1+s);
Lm(isnan(Lm)) = 0;
coef = Ls.*Lm;
end
function xk=idftc(cm)
cm=cm(:);
xk=cm;
km=length(cm)-1;
n=2*km+1;
mm=1:km;
kk=mm';
bkm=2*cos(2*pi/n*kk*mm);
cm2=cm(2:km+1);
xk(1)=(cm(1)+2*sum(cm2));
xk2=cm(1)+bkm*cm2;
xk=[xk(1);xk2];
end
