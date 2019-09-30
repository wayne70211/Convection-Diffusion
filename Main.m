%% Convection-Diffusion
clear
Max_time=8;
X_boundary=1;
dt=0.05;
dx=0.01;
X=0:dx:X_boundary;
t=0:dt:Max_time;
Nt=Max_time/dt+1;
Nx=X_boundary/dx+1;
u=0.08;
a=0.001;
CFL=u*dt/dx;
CF2L=a*dt/(dx*dx);
T_exact(1:Nt,1:Nx)=0;
tt=(Nt-1)/Max_time;

% Initial
T(1:Nt,1:Nx)=0;                                 % Create Space
T_n(1:Nx)=0;
XX=0:dx:0.2;
Nb=0.2/dx+1;
T(1,1:Nb)=1-(10.*(XX(1:Nb))-1).^2;              % Initial Condition

%% Exact Solution
for n=1:Nt
    for j=1:Nx
        if X(j)-u*t(n) <= 0.2 && X(j)-u*t(n) >= 0
            T_exact(n,j)=1-(10*(X(j)-u*t(n))-1)^2;
        else
            T_exact(n,j)=0;
        end
    end
end
 
%% Explicit Euler time advancement and second-order central difference for the spatial derivative.
for j=2:Nx-1
    T(2,j)=T(1,j)*(1-2*CF2L)+T(1,j-1)*(CFL/2+CF2L)+T(1,j+1)*(-CFL/2+CF2L);
end
 
for n=2:Nt-1
    for j=2:Nx-1
         T(n+1,j)=T(n,j)*(1-2*CF2L)+T(n,j-1)*(CFL/2+CF2L)+T(n,j+1)*(-CFL/2+CF2L);      
 %       T(n+1,j)=T(n,j)-CFL*(T(n,j+1)-T(n,j-1))/2+CF2L*(T(n-1,j+1)-2*T(n-1,j)+T(n-1,j-1)); %  Leapfrog on diffusion   
    end
end
 
figure(2)
subplot(3,1,1)
plot(X,T_exact(1,1:Nx),'linewidth',2)
title(['Convection-diffusion at ',num2str(0),' sec'])
hold on
grid on
plot(X,T(1,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution only Convection','Explicit Euler with Diffusion')
hold off
subplot(3,1,2)
plot(X,T_exact(4*tt,1:Nx),'linewidth',2)
title(['Convection-diffusion at ',num2str(4),' sec'])
hold on
grid on
plot(X,T(4*tt,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution only Convection','Explicit Euler with Diffusion')
hold off
subplot(3,1,3)
plot(X,T_exact(8*tt,1:Nx),'linewidth',2)
title(['Convection-diffusion at ',num2str(8),' sec'])
hold on
grid on
plot(X,T(8*tt,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution only Convection','Explicit Euler with Diffusion')

%% Leapfrog time advancement and the second-order central difference for the spatial derivative.
 
% Explicit Euler for first step
for j=2:Nx-1
    T(2,j)=T(1,j)-u*dt*(T(1,j+1)-T(1,j-1))/(2*dx);
end
 
for n=2:Nt-1
    for j=2:Nx-1
        T(n+1,j)=T(n-1,j)-u*dt*(T(n,j+1)-T(n,j-1))/(dx);
    end
end
 
%
figure(3)
subplot(3,1,1)
plot(X,T_exact(1,1:Nx),'linewidth',2)
title(['Pure convection at ',num2str(0),' sec'])
hold on
grid on
plot(X,T(1,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution','Leapfrog Solution')
hold off
subplot(3,1,2)
plot(X,T_exact(4*tt,1:Nx),'linewidth',2)
title(['Pure convection at ',num2str(4),' sec'])
hold on
grid on
plot(X,T(4*tt,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution','Leapfrog Solution')
hold off
subplot(3,1,3)
plot(X,T_exact(8*tt,1:Nx),'linewidth',2)
title(['Pure convection at ',num2str(8),' sec'])
hold on
grid on
plot(X,T(8*tt,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution','Leapfrog Solution')
 
%% Lax-Wendroff scheme
 
for n=1:Nt-1
    for j=2:Nx-1
        T(n+1,j)=T(n,j)-u*dt*(T(n,j+1)-T(n,j-1))/(2*dx)+(u*u)*dt*dt*(T(n,j+1)-2*T(n,j)+T(n,j-1))/(2*dx*dx);
    end
end
 
%
figure(4)
subplot(3,1,1)
plot(X,T_exact(1,1:Nx),'linewidth',2)
title(['Pure convection at ',num2str(0),' sec'])
hold on
grid on
plot(X,T(1,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution','Lax-Wendroff Solution')
hold off
subplot(3,1,2)
plot(X,T_exact(4*tt,1:Nx),'linewidth',2)
title(['Pure convection at ',num2str(4),' sec'])
hold on
grid on
plot(X,T(4*tt,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution','Lax-Wendroff Solution')
hold off
subplot(3,1,3)
plot(X,T_exact(8*tt,1:Nx),'linewidth',2)
title(['Pure convection at ',num2str(8),' sec'])
hold on
grid on
plot(X,T(8*tt,1:Nx),'linewidth',2)
xlabel('X')
ylabel('T')
legend('Exact Solution','Lax-Wendroff Solution')

%% Movie Maker
hfig=figure;
ax = gca;
axis tight
xlabel('X')
ylabel('T')
ax.NextPlot = 'replaceChildren';
v = VideoWriter('Simulation.avi','Motion JPEG AVI');
v.Quality = 95;
open(v);
for i=1:Nt-1
plot(X,T_exact(i,1:Nx),'linewidth',2)
hold on
plot(X,T(i,1:Nx),'linewidth',2)
title(['Solution at ',num2str(i*dt),' sec'])
axis([0,1,0,1.5])
grid on
hold off
drawnow
F(i) = getframe(hfig);
writeVideo(v,F(i))
end
close(v)
%% replay
figure(3)
movie(gcf,F)
