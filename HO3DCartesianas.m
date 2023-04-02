function M = HO3DCartesianas(nx,ny,nz,x,y,z,s)


[X,Y,Z]=meshgrid(x,y,z);

HGX=(2/(pi*s^2))^(1/4)*exp(-X.^2/s^2).*(1/sqrt(2^nx*factorial(nx))).*Hermite_pol(nx,sqrt(2).*X/s);
HGY=(2/(pi*s^2))^(1/4)*exp(-Y.^2/s^2).*(1/sqrt(2^ny*factorial(ny))).*Hermite_pol(ny,sqrt(2).*Y/s);
HGZ=(2/(pi*s^2))^(1/4)*exp(-Z.^2/s^2).*(1/sqrt(2^nz*factorial(nz))).*Hermite_pol(nz,sqrt(2).*Z/s);

M=HGX.*HGY.*HGZ;

end