function [coef] = ker(m,s)
Ls = gamma(s+0.5).*4.^s./(sqrt(pi).*abs(gamma(-s)));
Lm = gamma(abs(m)-s)./gamma(abs(m)+1+s);
Lm(isnan(Lm)) = 0;
coef = Ls.*Lm;
end