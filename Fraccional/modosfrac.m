function [normalfrac,V] = modosfrac(C,w0,n,cc,s)
% Devuelve las frecuencias y la matriz de autovectores V (los modos
% normales). Revisar si la matriz es la adecuada.
% cc = 1 Dirichlet
if cc == 1
   mat = zeros(n);
   for i = 2:n
   u = ones(1,n-i+1);
   mat = mat - diag(ker(i-1,s)*u,i-1);
   end
   mat = mat + mat';
   for i = 1:n
       mat(i,i) = -sum(mat(i,:));
   end
   mat = C*mat + w0^2*eye(n);
   [V,D] = eig(mat);
   normalfrac = sqrt(diag(D)); 
end
% Completar otras cc
end