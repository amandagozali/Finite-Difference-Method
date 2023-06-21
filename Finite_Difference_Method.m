clc;
clear all;
close all;

%Batas awal dan akhir
a = 0; 
b = pi;
t0 = 0;
ta = 5;

%Jumlah partisi
nx = 50;
nt = 750;

x = linspace(a,b,nx);
t = linspace(t0,ta,nt);

%Besar partisi
dx = (b-a)/nx;
dt = (ta-t0)/nt;

c = 1;
s = ((c.^2)*(dt.^2))/(dx.^2);


%Solusi analitik
%Inisialisasi matriks
u = zeros(nx,nt);

%Syarat awal
u(:,1) = sin(3.*x);

%Syarat Batas
u(1,:) = 0;
u(nx,:) = 0;

for i=1:nx
    for j = 1:nt
        u(i,j) = sin(3.*x(i))*cos(3.*t(j));
    end
end


%Solusi Beda Hingga Eksplisit
%Inisialisasi matriks
u_eks = zeros(nx,nt);

%Syarat awal
u_eks(:,1) = sin(3.*x);

%Syarat Batas
u_eks(1,:) = 0;
u_eks(nx,:) = 0;

%Perhitungan Solusi Eksplisit
for j = 1:nt-1
    for i = 2:nx-1 
        if j == 1
            u_eks(i,j+1) = s*(u_eks(i+1,j)-2*u_eks(i,j)+u_eks(i-1,j))+2*u_eks(i,j)-u_eks(i,j);
        else
            u_eks(i,j+1) = s*(u_eks(i+1,j)-2*u_eks(i, j)+u_eks(i-1,j))+2*u_eks(i,j)-u_eks(i,j-1);
        end
    end
end


%Membuat plot solusi eksplisit dan solusi analitik
for a  = 1:nt
    figure(1)
    plot(x,u(:,a),'b',x,u_eks(:,a),'r')
    title('Grafik Solusi Analitik dan Solusi Eksplisit')
    legend('Solusi Analitik','Solusi Eksplisit')
    xlim([0 pi])
    ylim([-1.5 1.5])
    pause(0.1)
end

t01 = 0*dt;
t02 = 374*dt;
t03 = 549*dt;
t04 = 749*dt;

figure;
plot(x,u(:,1),'b',x,u_eks(:,1),'r')
xlim([0 pi])
title(['Solusi pada t=',num2str(t01)])
legend('Solusi Analitik','Solusi Eksplisit','Location','north')
xlabel('x')
ylabel('u')

figure;
plot(x,u(:,375),'b',x,u_eks(:,375),'r')
xlim([0 pi])
title(['Solusi pada t=',num2str(t02)])
legend('Solusi Analitik','Solusi Eksplisit','Location','north')
xlabel('x')
ylabel('u')

figure;
plot(x,u(:,550),'b',x,u_eks(:,550),'r')
xlim([0 pi])
title(['Solusi pada t=',num2str(t03)])
legend('Solusi Analitik','Solusi Eksplisit','Location','north')
xlabel('x')
ylabel('u')

figure;
plot(x,u(:,750),'b',x,u_eks(:,750),'r')
xlim([0 pi])
title(['Solusi pada t=',num2str(t04)])
legend('Solusi Analitik','Solusi Eksplisit')
xlabel('x')
ylabel('u')

%Mencari galat
for aa=1:nt
    error = abs(u(:,aa)-u_eks(:,aa))./abs(u(:,aa));
    error(isnan(error)) = 0;
    er(aa) = sum(error)/nx;
end

%Membuat plot galat
figure;
plot(er)
title('Perhitungan Galat Metode Eksplisit')
xlabel('Partisi ke-(nt)')
ylabel('Galat')
text(25,15,'(79,4.15317)')
text(160,30,'(236,14.5321)')
text(330,60,'(393,40.1874)')
text(480,220,'(550,211.162)')
text(630,165,'(707,145.125)')

t1 = 78*dt;
t2 = 235*dt;
t3 = 392*dt;
t4 = 549*dt;
t5 = 706*dt;

figure;
subplot(2,3,1)
plot(x,u(:,79),'b',x,u_eks(:,79),'r')
title(['Solusi pada saat t=',num2str(t1)])
xlabel('x')
ylabel('u')

subplot(2,3,2)
plot(x,u(:,236),'b',x,u_eks(:,236),'r')
title(['Solusi pada saat t=',num2str(t2)])
xlabel('x')
ylabel('u')

subplot(2,3,3)
plot(x,u(:,393),'b',x,u_eks(:,393),'r')
title(['Solusi pada saat t=',num2str(t3)])
xlabel('x')
ylabel('u')

subplot(2,3,4.5)
plot(x,u(:,550),'b',x,u_eks(:,550),'r')
title(['Solusi pada saat t=',num2str(t4)])
xlabel('x')
ylabel('u')

subplot(2,3,5.5)
plot(x,u(:,707),'b',x,u_eks(:,707),'r')
title(['Solusi pada saat t=', num2str(t5)])
xlabel('x')
ylabel('u')


%Solusi Beda Hingga Implisit
%Syarat awal
U(:,1) = sin(3.*x);

%Syarat Batas
U(1,:) = 0;
U(nx,:) = 0;

%Inisialisasi matriks
m = zeros(nx,nx);

for i = 1:nx
    for k = 1:nx
        if k == i-1
            m(i,k) = -s;
        elseif k == i 
            m(i,k) = 1+2*s;
        elseif k == i+1
            m(i,k) = -s;
        end
    end
end

U1 = inv(m)*U;
A = [U U1];

for i = 1 : nt - 1
   U2 = inv(0.5*m)*(U1-(0.5*U));
   u_imp = [A U2];
   U = U1;
   U1 = U2;
   A = u_imp;
end


%Membuat plot solusi implisit dan solusi analitik
for b = 1:nt
    figure(8)
    plot(x,u(:,b),'b',x, u_imp(:, b),'g')
    title('Grafik Solusi Analitik dan Solusi Implisit')
    legend('Solusi Analitik','Solusi Implisit')
    xlim([0 pi])
    ylim([-1.5 1.5])
    pause(0.1)
end

t01 = 0*dt;
t02 = 374*dt;
t03 = 549*dt;
t04 = 749*dt;

figure;
plot(x,u(:,1),'b',x,u_imp(:,1),'g')
xlim([0 pi])
title(['Solusi pada t=',num2str(t01)])
legend('Solusi Analitik','Solusi Implisit','Location','north')
xlabel('x')
ylabel('u')

figure;
plot(x,u(:,375),'b',x,u_imp(:,375),'g')
xlim([0 pi])
title(['Solusi pada t=',num2str(t02)])
legend('Solusi Analitik','Solusi Implisit','Location','north')
xlabel('x')
ylabel('u')

figure;
plot(x,u(:,550),'b',x,u_imp(:,550),'g')
xlim([0 pi])
title(['Solusi pada t=',num2str(t03)])
legend('Solusi Analitik','Solusi Implisit')
xlabel('x')
ylabel('u')

figure;
plot(x,u(:,750),'b',x,u_imp(:,750),'g')
xlim([0 pi])
title(['Solusi pada t=',num2str(t04)])
legend('Solusi Analitik','Solusi Implisit')
xlabel('x')
ylabel('u')

%Mencari galat 
for bb=1:nt
    error = abs(u(:,bb)-u_imp(:,bb));
    err(bb) = sum(error)/nx;
end

%Membuat plot galat
figure;
plot(err)
title('Perhitungan Galat Metode Implisit')
xlabel('Partisi ke-(nt)')
ylabel('Galat')
xlim([0 750])
text(20,0.06,'(73,0.0429369)')
text(100,0.11,'(180,0.0933362)')
text(240,0.13,'(292,0.119172)')
text(380,0.15,'(425,0.13233)')
text(510,0.18,'(576,0.164109)')
text(600,0.22,'(743,0.211352)')

t11 = 72*dt;
t22 = 179*dt;
t33 = 291*dt;
t44 = 424*dt;
t55 = 575*dt;
t66 = 742*dt;

figure;
subplot(2,3,1)
plot(x,u(:,73),'b',x,u_imp(:,73),'g')
title(['Solusi pada saat t=',num2str(t11)])
xlabel('x')
ylabel('u')

subplot(2,3,2)
plot(x,u(:,180),'b',x,u_imp(:,180),'g')
title(['Solusi pada saat t=',num2str(t22)])
xlabel('x')
ylabel('u')

subplot(2,3,3)
plot(x,u(:,292),'b',x,u_imp(:,292),'g')
title(['Solusi pada saat t=',num2str(t33)])
xlabel('x')
ylabel('u')

subplot(2,3,4)
plot(x,u(:,425),'b',x,u_imp(:,425),'g')
title(['Solusi pada saat t=',num2str(t44)])
xlabel('x')
ylabel('u')

subplot(2,3,5)
plot(x,u(:,576),'b',x,u_imp(:,576),'g')
title(['Solusi pada saat t=',num2str(t55)])
xlabel('x')
ylabel('u')

subplot(2,3,6)
plot(x,u(:,743),'b',x,u_imp(:,743),'g')
title(['Solusi pada saat t=',num2str(t66)])
xlabel('x')
ylabel('u')