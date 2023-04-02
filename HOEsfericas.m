function M = HOEsfericas(n,l,m,x,y,z,s)

[X,Y,Z]=meshgrid(x,y,z);

r=sqrt(X.^2+Y.^2+Z.^2);
theta=atan2(sqrt(X.^2+Y.^2),Z);
psi=atan2(Y,X);

const=sqrt(2^(l+5/2)*factorial(n)/(s^3*gamma(n+l+3/2)));
Radio = const.*(r/s).^l.*exp(-r.^2/s.^2).*PolyLaguerre(n,l+(1/2),(2*r.^2/s^2));

Plm = legendre(l,cos(theta));
if l ~= 0
    Plm = reshape(Plm(abs(m)+1,:,:),size(psi));
end
Ylm=sqrt((2*l+1)/(4*pi)*factorial(l-abs(m))/factorial(l+abs(m)))*Plm.*exp(1i*m*psi);

M=Radio.*Ylm;

end