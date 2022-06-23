function [mat] = matfrac(n,cc,s)
% Este código construye la matriz del núcleo, no multiplicamos por la C dentro 
% El término diagonal es la suma del resto. Es válida para la ecuación
% diferencial.
% cc = 1, El extremo está libre (sólo nota fuerzas por un lado)
if cc == 1
   mat = zeros(n);
   for i = 2:n
   u = ones(1,n-i+1);
   mat = mat - diag(ker(i-1,s)*u,i-1);
   %+ diag(ker(i-1,s)*v);
   end
   mat = mat + mat';
   for i = 1:n
       mat(i,i) = -sum(mat(i,:)); %Esto es lo correcto en extremos libres
   end
   %normales = sqrt(eig(mat)); 
end
% Completar otras cc
end