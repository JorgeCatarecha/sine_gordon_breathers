function [normales,mat] = modos(C,w0,n,cc)
% cc = 1 Dirichlet
if cc == 1
   u = ones(1,n);
   v = ones(1,n-1);
   mat = diag((w0^2+2*C)*u)-diag(C*v,1)-diag(C*v,-1);
   normales = sqrt(eig(mat)); 
end
% cc = 0 Periodicas
if cc == 0
   u = ones(1,n);
   v = ones(1,n-1);
   mat = diag((w0^2+2*C)*u)-diag(C*v,1)-diag(C*v,-1);
   mat(1,n) = -C;
   mat(n,1) = -C;
   normales = sqrt(eig(mat)); 
end
end