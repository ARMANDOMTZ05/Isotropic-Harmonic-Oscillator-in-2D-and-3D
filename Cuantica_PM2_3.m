%% Problema 3.4
clc
clear
close all

r=linspace(0,5,200);

r1=Radio3DHO(2,1,r,1);

c_cons1=trapz(r,r1.*conj(r1).*r.^2);

r2=Radio3DHO(1,2,r,1);

c_cons2=trapz(r,r2.*conj(r2).*r.^2);

figure(1)
plot(r,r1)
hold on
plot(r,r2)
ylabel('$R(r)$',Interpreter='latex')
xlabel('$r$',Interpreter='latex')
title('Variacion radial',Interpreter='latex')
legend('\psi_{2,1}^1','\psi_{1,2}^0')

%% Problema 3.5
%Ejercico 2

x=linspace(-3,3);
y=linspace(-3,3);
z=linspace(-3,3);

M1=HO3DCartesianas(0,1,2,x,y,z,1);
M2=HOEsfericas(0,2,-1,x,y,z,1);

index_cart=[0 0 0;1 0 0; 0 1 0; 0 0 1;1 1 0; 0 1 1;1 0 1;2 0 0; 0 2 0; 0 0 2;1 1 1;2 1 0;0 2 1; 1 0 2;2 0 1;0 1 2;1 2 0;3 0 0;0 3 0; 0 0 3; ...
    2 1 1;1 2 1;1 1 2; 2 2 0; 0 2 2;2 0 2;3 1 0;3 0 1;0 3 1;1 3 0; 0 1 3;1 0 3;4 0 0;0 4 0; 0 0 4];
index_esf=[0 0 0; 0 1 -1; 0 1 0; 0 1 1;1 0 0;0 2 -2;0 2 -1; 0 2 0;0 2 1; 0 2 2;1 1 -1;1 1 0; 1 1 1;0 3 -3;0 3 -2; 0 3 -1; 0 3 0; 0 3 1; 0 3 2; 0 3 3;...
    2 0 0;1 2 -2;1 2 -1; 1 2 0; 1 2 1; 1 2 2; 0 4 -4; 0 4 -3; 0 4 -2; 0 4 -1; 0 4 0; 0 4 1; 0 4 2; 0 4 3;0 4 4];

aprox_cart=zeros(length(x),length(y),length(z));

%Aproximación de esférica con cartesianas
for ii=1:35
    psi=HO3DCartesianas(index_cart(ii,1),index_cart(ii,2),index_cart(ii,3),x,y,z,1);
    c_n(ii)=trapz(x,trapz(y,trapz(z,conj(psi).*M2,3),2));
    aprox_cart=aprox_cart+c_n(ii)*psi;
end

aprox_esf=zeros(length(x),length(y),length(z));

%Aproximacion de cartesiana con esféricas
for jj=1:35
    psi2=HOEsfericas(index_esf(jj,1),index_esf(jj,2),index_esf(jj,3),x,y,z,1);
    c_n2(jj)=trapz(x,trapz(y,trapz(z,conj(psi2).*M1,3),2));
    aprox_esf=aprox_esf+c_n2(jj)*psi2;
end




%% Problema 3.6

%Esféricas aproximada
serpiente1(:,:)=aprox_cart(:,50,:);
serpiente2(:,:)=aprox_cart(50,:,:);

%Cartesiana aproximada
serpiente3(:,:)=aprox_esf(:,50,:);
serpiente4(:,:)=aprox_esf(50,:,:);

%Esfereica real
serpiente5(:,:)=M2(:,50,:);
serpiente6(:,:)=M2(50,:,:);

%Cartesiana real
serpiente7(:,:)=M1(:,50,:);
serpiente8(:,:)=M1(50,:,:);

figure(1)
subplot(3,2,1)
imagesc(x,y,abs(aprox_esf(:,:,50)).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Aproximacion de $\psi_{0,1}^{2}$ con eigenestados esfericos en $z=0$',Interpreter='latex')
subplot(3,2,2)
imagesc(x,y,abs(M1(:,:,50)).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestado esferico $\psi_{0,1}^{2}$ en $z=0$',Interpreter='latex')
subplot(3,2,3)
imagesc(x,y,abs(serpiente3).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$z$',Interpreter='latex')
title('Aproximacion de $\psi_{0,1}^{2}$ con eigenestados esfericos en $y=0$',Interpreter='latex')
subplot(3,2,4)
imagesc(x,y,abs(serpiente7).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$z$',Interpreter='latex')
title('Eigenestado esferico $\psi_{0,1}^{2}$ en $y=0$',Interpreter='latex')
subplot(3,2,5)
imagesc(x,y,abs(serpiente4).^2)
colormap jet
axis equal
xlabel('$z$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Aproximacion de $\psi_{0,1}^{2}$ con eigenestados esfericos en $x=0$',Interpreter='latex')
subplot(3,2,6)
imagesc(x,y,abs(serpiente8).^2)
colormap jet
axis equal
xlabel('$z$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestado esferico $\psi_{0,1}^{2}$ en $x=0$',Interpreter='latex')

figure(2)
subplot(3,2,1)
imagesc(x,y,abs(aprox_cart(:,:,50)).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Aproximacion de $\psi_{0,2}^{-1}$ con eigenestados esfericos en $z=0$',Interpreter='latex')
subplot(3,2,2)
imagesc(x,y,abs(M2(:,:,50)).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestado $\psi_{0,2}^{-1}$ en $z=0$',Interpreter='latex')
subplot(3,2,3)
imagesc(x,y,abs(serpiente1).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$z$',Interpreter='latex')
title('Aproximacion de $\psi_{0,2}^{-1}$ con eigenestados esfericos en $y=0$',Interpreter='latex')
subplot(3,2,4)
imagesc(x,y,abs(serpiente5).^2)
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$z$',Interpreter='latex')
title('Eigenestado $\psi_{0,2}^{-1}$ en $y=0$',Interpreter='latex')
subplot(3,2,5)
imagesc(x,y,abs(serpiente2).^2)
colormap jet
axis equal
xlabel('$z$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Aproximacion de $\psi_{0,2}^{-1}$ con eigenestados esfericos en $x=0$',Interpreter='latex')
subplot(3,2,6)
imagesc(x,y,abs(serpiente6).^2)
colormap jet
axis equal
xlabel('$z$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestado $\psi_{0,2}^{-1}$ en $x=0$',Interpreter='latex')




%% 
clc
clear
close all

x=linspace(-3,3,200);
y=linspace(-3,3,200);
z=linspace(-3,3,200);

psi=(1/3)*HOEsfericas(0,2,1,x,y,z,1)-(2i/5)*HOEsfericas(1,1,-1,x,y,z,1)+0.5*HOEsfericas(1,2,2,x,y,z,1);
c_cost=sqrt(1/trapz(x,trapz(y,trapz(z,psi.*conj(psi),3),2)));

norm_psi=c_cost*psi;

I=trapz(x,trapz(y,trapz(z,norm_psi.*conj(norm_psi),3),2));

serpiente(:,:)=norm_psi(:,100,:);
serpiente2(:,:)=norm_psi(100,:,:);

figure(2)
subplot(2,2,1)
imagesc(x,y,abs(norm_psi(:,:,100)).^2);
colormap jet
colorbar
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Plano $XY$ en $z=0$',Interpreter='latex')

subplot(2,2,2)
imagesc(x,z,abs(serpiente).^2)
colormap jet
xlabel('$x$',Interpreter='latex')
ylabel('$z$', Interpreter='latex')
title('Plano $XZ$ en $y=0$',Interpreter='latex')
colorbar

subplot(2,2,[3,4])
imagesc(y,z,abs(serpiente2).^2)
colormap jet
xlabel('$y$',Interpreter='latex')
ylabel('$z$',Interpreter='latex')
title('Plano $YZ$ en $x=0$',Interpreter='latex')
colorbar


M1=(1/3)*SphArmonic(2,1)-(2i/5)*SphArmonic(1,-1)+0.5*SphArmonic(2,2);

az=linspace(0,2*pi,200);
col=linspace(0,pi,200);


[phi,theta] = meshgrid(az,col);

I2=trapz(az,trapz(col,abs(M1).^2.*sin(theta),2));

[xn,yn,zn]=sphere(199);

phi2=phi-pi;

x_hammer=(sqrt(8).*sin(theta).*sin(phi2/2))./sqrt(1+sin(theta).*cos(phi2/2));
y_hammer=(sqrt(2).*cos(theta))./sqrt(1+sin(theta).*cos(phi2/2));

figure(3)
subplot(2,1,1)
surf(xn,yn,zn,abs(M1).^2.*sin(theta),'EdgeColor','none')
axis equal
colormap jet
colorbar
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
zlabel('$z$',Interpreter='latex')
title('Funcion de densidad de probabilidad $|\Phi|^2$', Interpreter='latex')
subplot(2,1,2)
surf(x_hammer,y_hammer,abs(M1).^2,'EdgeColor','none')
xlabel('x')
ylabel('y')
colormap jet
colorbar
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
zlabel('$z$',Interpreter='latex')
title('Funcion de densidad de probabilidad $|\Phi|^2$ en proyeccion de Hammer', Interpreter='latex')
view(2)
axis equal
