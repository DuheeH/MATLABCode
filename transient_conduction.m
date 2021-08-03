function transient_conduction
%Created by Miguel Modestino 10/01/2018
%Example of Transient Conduction solved using ODE solver
%Chemical and Biomolecular Engineering
%New York University
%Heat and Mass Transport


Ti = 25; % C

L = 0.01; % m

N = 20; % Number of grid points

T0 = ones(1,N)*Ti; %Initial Condition

[t_out,T_out]=ode15s(@ode_transient,[0 1000],T0); %ODE solver

x = linspace(0,L,N); 

surf(x,t_out,T_out)
xlabel('x [m]')
ylabel('time [s]')
zlabel('Temperature [C]')

end


function dTdt=ode_transient(t,T)
%ODEFUN

k = 0.1;% W/m-K
rho = 8500;%Kg/m3
c = 400; %400 J/kg-K
h = 400; %W/m^2-K

alpha = k/rho/c; %thermal diffusivity

Tinf = 200; % C
Ti = 25; %25 C

L = 0.01; % m

N = 20; % Number of grid points

dx = L/(N-1);

dTdt=zeros(N,1);
%First B.C.
dTdt(1) = alpha/dx^2*(-T(1)+T(2));

%Internal Points
for j = 2:N-1
    dTdt(j)=alpha/dx^2*(T(j-1)-2*T(j)+T(j+1));
end

%Final Boundary Condition
%dTdt(N)=alpha/dx^2*(T(N-1)-2*T(N) + 1/(h+k/dx)*(h*Tinf-k/dx*T(N)));
dTdt(N)=alpha/dx^2*(T(N-1) -(1+h*dx/k)*T(N)+ h*dx/k*Tinf);

end
