%% Proyecto Mensaual 2

%Problema 1.1
clc
clear
close all

r=0:0.1:30;
a=1;
[Prob1,int1]=FuncAtomoHidrogeno(r,a,1,0);
[Prob2,int2]=FuncAtomoHidrogeno(r,a,2,0);
[Prob3,int3]=FuncAtomoHidrogeno(r,a,2,1);
[Prob4,int4]=FuncAtomoHidrogeno(r,a,3,0);
[Prob5,int5]=FuncAtomoHidrogeno(r,a,3,1);
[Prob6,int6]=FuncAtomoHidrogeno(r,a,3,2);

figure(1)
plot(r,Prob1,r,Prob2,r,Prob3,r,Prob4,r,Prob5,r,Prob6)
title('Probabilidad de la funcion radial de atomo de Hidrogeno',Interpreter='latex')
xlabel('$r$',Interpreter='latex')
ylabel('$r^2|R^l_n(r)|^2$',Interpreter='latex')
legend('n=1, l=0','n=2, l=0','n=2, l=1','n=3, l=0','n=3, l=1','n=3, l=2')


%% Parte 3
%Problema 2.1.1

clc
clear
close all

x=-3:0.01:3;
y=-3:0.01:3;
n_x=1;
n_y=2;

h=1;
M=2;
w=1;
s=sqrt(2*h/(M*w));
phi=HO_2Dimension(n_x,n_y,x,y,s);
phi2=HO_2Dimension(3,2,x,y,s);
figure(3)
subplot(2,1,1)
surf(x,y,phi.*phi,'EdgeColor','none');
colorbar
colormap jet
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestados Cartesianos $\psi_{1,2}$',Interpreter='latex')
view(2)

subplot(2,1,2)
surf(x,y,phi2.*phi2,'EdgeColor','none');
colorbar
colormap jet
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestados Cartesianos $\psi_{3,2}$',Interpreter='latex')
view(2)


%% Oscilador armÛnico en dos dimensiones Estados cartesianos
%Problema 2.1,2 al 4

clc
clear
close all

x=-5:0.1:5;
y=-5:0.1:5;



h=1;
M=2;
w=1;
s=sqrt(2*h/(M*w));
En=[2*h*w;5*h*w;4*h*w];

phi_1=HO_2Dimension(0,1,x,y,s);
phi_2=HO_2Dimension(3,1,x,y,s);
phi_3=HO_2Dimension(1,2,x,y,s);

psi=(1/3)*phi_1-(2i/5)*phi_2+(1/2)*phi_3;


c=sqrt(1/trapz(x,trapz(y,abs(psi).^2,2)));

norm_psi=c*psi;


figure(5)
surf(x,y,abs(norm_psi).^2,EdgeColor="none")
colormap jet
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
zlabel('$|\Psi|^2$',Interpreter='latex')
title('$C\left(\frac{1}{2}\psi_{0,1}-\frac{i2}{5}\psi_{3}{1}+\frac{1}{2}\psi{1}{2}\right)$',Interpreter='latex')
view(2)
axis equal
colorbar



arg=angle(norm_psi);

figure(6)
surf(x,y,arg,EdgeColor="none")
colormap jet
title('Fase de $C\left(\frac{1}{2}\psi_{0,1}-\frac{i2}{5}\psi_{3}{1}+\frac{1}{2}\psi{1}{2}\right)$ ','Interpreter','latex')
view(2)
axis equal
colorbar

[X,Y]=meshgrid(x,y);
[GX,GY]=gradient(norm_psi);
Jx=(h/M)*imag(conj(norm_psi).*GX);
Jy=(h/M)*imag(conj(norm_psi).*GY);



figure(7)
quiver(X,Y,Jx,Jy)
axis([-2 2 -2 2])
title('Flujo $\frac{\hbar}{M}Im(\psi^*\nabla\psi)$','Interpreter','latex')
colorbar


time=linspace(0,2*pi/w);
t=zeros(1,1,length(time));
t(1,1,:)=time(:);

% for k=1:length(t)
%     psi_time(:,:,k)=c*((1/3)*phi_1*exp(1i*En(1)*t(1,k))-(2i/5)*phi_2*exp(1i*En(2)*t(1,k))+(1/2)*phi_3*exp(1i*En(3)*t(1,k)));
% end

psi_time=c*((1/3)*phi_1.*exp(1i*En(1).*t)-(2i/5)*phi_2.*exp(1i*En(2).*t)+(1/2)*phi_3.*exp(1i*En(3).*t));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
abs_psi_time=psi_time.*conj(psi_time);

figure(8)
for ii=1:length(time)
    surf(x,y,abs_psi_time(:,:,ii),EdgeColor="none")
    colormap jet
    xlabel('$x$',Interpreter='latex')
    ylabel('$y$',Interpreter='latex')
    zlabel('$|\Psi|^2$','Interpreter','latex')
    title('$\frac{1}{3}\psi_{01}-\frac{2i}{5}\psi_{31}+\psi_{12}$',Interpreter='latex')
    axis([-5 5 -5 5 0 1]);

    pause(0.05)
end
%Animar esta parte

%% Oscilador armÛnico en dos dimensiones Circulares

 %Problema 2.2.1

clc
clear
close all

x=linspace(-5,5);
y=linspace(-5,5);

h=1;
M=2;
w=1;
s=sqrt(2*h/(M*w));

phi_circular=HO_circular(2,1,x,y,s);
phi_circular2=HO_circular(3,2,x,y,s);

orto=real(phi_circular.*phi_circular2);

dx=trapz(x,orto,1);
dy=trapz(y,dx,2);

arg2=angle(phi_circular);
arg3=angle(phi_circular2);

figure(8)
subplot(2,2,1)
surf(x,y,phi_circular.*conj(phi_circular),EdgeColor='none');
colormap jet
axis equal
view(2)
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestados Polares $\psi_{2,1}$',Interpreter='latex')

subplot(2,2,2)
surf(x,y,phi_circular2.*conj(phi_circular2),EdgeColor='none')
colormap jet
axis equal
view(2)
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestados Polares $\psi_{3,2}$',Interpreter='latex')


subplot(2,2,3)
surf(x,y,arg2,EdgeColor="none")
colormap jet
axis equal
view(2)
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Fase Eigenestados Polares $\psi_{2,1}$',Interpreter='latex')

subplot(2,2,4)
surf(x,y,arg3,EdgeColor="none")
colormap jet
axis equal
view(2)
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Fase Eigenestados Polares $\psi_{3,2}$',Interpreter='latex')


%% Problema 2.2.2

clc
clear
close all

x=linspace(-5,5,200);
y=linspace(-5,5,200);

h=1;
M=2;
w=1;
s=sqrt(2*h/(M*w));

psi_circular=HO_circular(0,2,x,y,s);
psi_circular2=HO_circular(1,-1,x,y,s);
psi_circular3=HO_circular(2,0,x,y,s);


psi=(1/3)*psi_circular-(2i/5)*psi_circular2+(1/2)*psi_circular3;
% 
% figure(12)
% graph10=surf(x,y,psi.*conj(psi));
% graph10.EdgeColor='none';
% 
Q=trapz(x,psi.*conj(psi),1);
Q2=trapz(y,Q,2);
c=sqrt(1/Q2);

norm_psi=c*psi;

I=trapz(y,trapz(x,norm_psi.*conj(norm_psi),2));

arg=angle(norm_psi);

figure(13)
graph11=surf(x,y,norm_psi.*conj(norm_psi));
graph11.EdgeColor='none';
colormap jet
colorbar
axis equal
view(2)
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('$|\Phi\rangle=C\left(\frac{1}{3}\psi_0^2-\frac{2i}{5}\psi_1^{-1}+\frac{1}{2}\psi_2^0\right)$','Interpreter','latex')
colorbar

figure(14)
surf(x,y,arg,EdgeColor="none")
colormap jet
axis equal
view(2)
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Fase $|\Phi\rangle=C\left(\frac{1}{3}\psi_0^2-\frac{2i}{5}\psi_1^{-1}+\frac{1}{2}\psi_2^0\right)$',Interpreter='latex')
colorbar

%% Problema 2.3.1
clc
clear
close all

tic

x=linspace(-5,5);
y=linspace(-5,5);

n_x=0:1:5;
n_y=0:1:5;
c_n=zeros(length(n_x),length(n_y));

h=1;
M=2;
w=1;
s=sqrt(2*h/(M*w));

psi_circular=HO_circular(2,1,x,y,s);

TD_gauss_eq=zeros(length(x),length(y));


for k=1:length(n_x)
    for l=1:length(n_y)
        psi_n=HO_2Dimension(n_x(k),n_y(l),x,y,s);
        c_n(k,l)=trapz(y,trapz(x,conj(psi_n).*psi_circular,2));
        TD_gauss_eq=TD_gauss_eq+c_n(k,l)*psi_n;

    end
end



c_const=sqrt(1/trapz(y,trapz(x,conj(TD_gauss_eq).*TD_gauss_eq,2)));
norm_approx_cart=c_const*TD_gauss_eq;

%% 

figure (14)
subplot(2,1,1)
surf(x,y,psi_circular.*conj(psi_circular),'EdgeColor','none')
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Eigenestado Polar $\psi_2^1$',Interpreter='latex')
view(2)
colorbar

subplot(2,1,2)
surf(x,y,conj(norm_approx_cart).*norm_approx_cart,'EdgeColor','none')
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Aproximacion con eigenestados cartesianos',Interpreter='latex')
view(2)
colorbar
toc

%% Pulso Gaussiano Problema 2.3.2
clc
clear
close all

x=linspace(-5,5,200);
y=linspace(-5,5,200);

h=1;
M=2;
w=1;
s=sqrt(2*h/(M*w));

[X,Y]=meshgrid(x,y);

psi_gauss=exp(-(X-s).^2/s^2).*exp(-(Y+1.5*s).^2/(3*s/5)^2).*exp(2i*pi/s*(X/3+Y/4));


c=sqrt(1/trapz(x,trapz(y,conj(psi_gauss).*psi_gauss,2)));

norm_psi_gauss=c*psi_gauss;

integral=trapz(y,trapz(x,conj(norm_psi_gauss).*norm_psi_gauss,2));


[GX,GY]=gradient(norm_psi_gauss);

psi_gauss_momentum=real(conj(norm_psi_gauss).*(-1i*h*(GX+GY)));

Momentum=trapz(y,trapz(x,psi_gauss_momentum,2));

psi_gauss_AMomentum=real(conj(norm_psi_gauss).*(-1i*h*(X.*GY-Y.*GX)));

AMomentum=trapz(y,trapz(x,psi_gauss_AMomentum,2));


n_x=0:1:10;
n_y=0:1:10;
c_n=zeros(length(n_x),length(n_y));

% phi_gauss_x=exp(-(x-s).^2/s^2).*exp(2i*pi/s*(x/3));
% phi_gauss_y=exp(-(y+1.5*s).^2/(3*s/5)^2).*exp(2i*pi/s*(y/4));

E_n=(n_x'+n_y+1)*h*w;

time=linspace(0,2*pi/w);

TD_guss_eq=zeros(length(x),length(y),length(time));

TD_gauss_eq=zeros(length(x),length(y));

t(1,1,:)=time(:);

% for k=1:length(n_x)
%     for l=1:length(n_y)
%         [~,~,psi_n]=HO_2Dimension(n_x(k),n_y(l),x',y,s);
%         c_n(k,l)=trapz(y,trapz(conj(psi_n).*norm_psi_gauss,2));
%         approx_carterianas=approx_carterianas+c_n(k,l)*psi_n.*exp(-1i*E_n(k,l).*t/h);
%     end
% 
% end
% 
% c_const=sqrt(1/trapz(y,trapz(x,conj(approx_carterianas(:,:,1)).*approx_carterianas(:,:,1),2)));
% norm_approx_cart=c_const*approx_carterianas;
% 
% 
% figure(17)
% surf(x,y,conj(norm_approx_cart(:,:,1)).*norm_approx_cart(:,:,1),'EdgeColor','none')


m=-10:1:10;
n=0:1:20;
E_q=(2*n'+abs(m)+1)*h*w;

for k=1:length(m)
    for l=1:length(n)
        psi_n=HO_circular(n(k),m(l),x,y,s);
        c_n(k,l)=trapz(y,trapz(conj(psi_n).*norm_psi_gauss,2));
        TD_gauss_eq=TD_gauss_eq+c_n(k,l)*psi_n.*exp(-1i*E_q(k,l).*t/h);
    end
end

circular_c=sqrt(1/trapz(y,trapz(x,conj(TD_gauss_eq(:,:,1)).*TD_gauss_eq(:,:,1),2)));
norm_approx_circular=circular_c*TD_gauss_eq;

[X,Y]=meshgrid(x,y);
x_expected=trapz(y,trapz(x,conj(norm_approx_circular).*norm_approx_circular.*X,2));
rx=zeros(1,length(x_expected));
rx(1,:,1)=x_expected;
y_expected=trapz(y,trapz(x,conj(norm_approx_circular).*norm_approx_circular.*Y,2));
ry=zeros(1,length(y_expected));
ry(1,:,1)=y_expected;

%% 


figure(18)
subplot(2,1,1)
surf(x,y,abs(norm_approx_circular(:,:,1)).^2,EdgeColor="none")
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Aproximacion con eigenestados circulares',Interpreter='latex')
view(2)
colorbar

subplot(2,1,2)
surf(x,y,abs(norm_psi_gauss).^2,'EdgeColor','none')
colormap jet
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Pulso Gaussiano',Interpreter='latex')
view(2)
colorbar



figure(2)
plot(rx,ry,'r')
axis equal
xlabel('$x$',Interpreter='latex')
ylabel('$y$',Interpreter='latex')
title('Valor esperado de la posicion $\langle \vec{r}\rangle$',Interpreter='latex')
view(2)


figure(3)
B = bar3(m,abs(circular_c*c_n));
xlabel('$n$',Interpreter='latex')
ylabel('$m$',Interpreter='latex')
title('$|C_{n,m}|$',Interpreter='latex')
colorbar
view(2)
colormap jet
for k = 1:length(B)
     zdata = B(k).ZData;
     B(k).CData = zdata;
     B(k).FaceColor = 'interp';
end
iter=length(n);
ticks=cell(iter,1);
for tick=1:iter
    ticks{tick}=string(tick-1);
end

xticklabels(7.*ticks)


% figure(19)
% for ii = 1:length(time)
% 
%     surf(x,y,conj(norm_approx_circular(:,:,ii)).*norm_approx_circular(:,:,ii),EdgeColor="none")
%     colormap jet
%     hold on
%     plot(rx(1:ii),ry(1:ii),'r.','LineWidth', 2)
%     hold off
%     xlabel('$x$',Interpreter='latex')
%     ylabel('$x$',Interpreter='latex')
%     zlabel('|\Psi|^2')
%     title('Aproximacion con eigenestados polares $\psi_n^l$',Interpreter='latex')
%     axis([-5 5 -5 5 0 1]);
%     view(2)
% 
%     pause(0.05)
% end


% figure(20)
% for ii = 1:length(time)
% 
%     surf(x,y,conj(norm_approx_cart(:,:,ii)).*norm_approx_cart(:,:,ii),EdgeColor="none")
%     xlabel('x')
%     ylabel('y')
%     zlabel('|\Psi|^2')
%     axis([-5 5 -5 5 0 1]);
% 
%     pause(0.05)
% end
