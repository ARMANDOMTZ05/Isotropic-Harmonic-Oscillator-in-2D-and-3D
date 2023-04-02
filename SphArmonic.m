function Armonico = SphArmonic(l,m)

% [X,Y,Z]=meshgrid(x,y,z);
% 
% const=sqrt((2*l+1)/(4*pi)*factorial(l-abs(m))/factorial(l+abs(m)));
% 
% theta=atan2(Y,X);
% 
% psi=atan2(sqrt(X.^2+Y.^2),Z);
% 
% Plm = legendre(l,cos(theta));
% if l ~= 0
%     Plm = reshape(Plm(m+1,:,:),size(X));
% end
% 
% Armonico=const*Plm.*exp(1i*m*psi);
% 

az=linspace(0,2*pi,200);
col=linspace(0,pi,200);

[phi,theta] = meshgrid(az,col);

Plm = legendre(l,cos(theta));

if l ~= 0
    Plm = reshape(Plm(abs(m)+1,:,:),size(phi));
end

a = (2*l+1)*factorial(l-abs(m));
b = 4*pi*factorial(l+abs(m));
C = sqrt(a/b);
Armonico = C .*Plm .*exp(1i*m*phi);

% = sph2cart(phi, pi/2-theta, abs(real(Armonico)));

end