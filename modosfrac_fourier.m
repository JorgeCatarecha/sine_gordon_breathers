function [normalfrac] = modosfrac_fourier(C,w0,n,s)

m=toeplitz([0:n/2-1 n/2:-1:1]); %Problema con condiciones de contorno periódicas

% Aquí se introduce el Kernel K(n)

% Fraccional
Kernel = ker(m,s);
Kernel = Kernel - diag(diag(Kernel)); %Para quitar el término central
Kernel(isinf(Kernel))=0; %Por si el kernel diverge
mat=Kernel-sum(Kernel(1,:))*eye(n);

   mat = -C*mat + w0^2*eye(n);
   [~,D] = eig(mat);
   normalfrac = sqrt(diag(D)); 
   function [coef] = ker(m,s)
        Ls = gamma(s+0.5)*4^s./(sqrt(pi)*abs(gamma(-s)));
        Lm = gamma(abs(m)-s)./gamma(abs(m)+1+s);
        Lm(isnan(Lm)) = 0;
        coef = Ls.*Lm;
    end
end