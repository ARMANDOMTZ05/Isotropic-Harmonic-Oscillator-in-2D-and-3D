function Radio = Radio3DHO(n,l,r,s)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
const=sqrt(2^(l+5/2)*factorial(n)/(s^3*gamma(n+l+3/2)));
Radio = const.*(r/s).^l.*exp(-r.^2/s.^2).*PolyLaguerre(n,l+(1/2),(2*r.^2/s^2));
end