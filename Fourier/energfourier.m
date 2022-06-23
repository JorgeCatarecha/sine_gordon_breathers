function [E] = energfourier(z,C,alfa)
% Calcula la energía de la solución u, con una constante de acoplo C
% cc periodicas
[coef_lim,~] = size(z);
pos=ones(1,coef_lim)*[z(1,:);2*z(2:end,:)];
n = length(pos);
m=ker(toeplitz([0:n/2-1 n/2:-1:1]),alfa);
U1 = 0;
for i=1:n
U1 = U1 + energiapot(i); %Energía potencial elástica
end
U2 = sum(1-cos(pos)); %Energía potencial tipo coseno
E = C*U1/2 + U2;
    function [energia] = energiapot(i)
        coef = m(i,:);
        energia = coef*(pos(i)-pos').^2;
    end
function [coef] = ker(m,s)
        Ls = gamma(s+0.5)*4^s./(sqrt(pi)*abs(gamma(-s)));
        Lm = gamma(abs(m)-s)./gamma(abs(m)+1+s);
        Lm(isnan(Lm)) = 0;
        coef = Ls.*Lm;
    end
end