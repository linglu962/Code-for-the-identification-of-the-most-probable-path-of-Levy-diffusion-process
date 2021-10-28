clear;clc;
global f h Mf I diagg zz A B C D;
% mesh parameters
zmin=220;zmax=320;
I=1000;h=1/I;
x=linspace(-1+h,1-h,2*I-1);z=(zmax-zmin)/2*x+(zmax+zmin)/2;nx=size(x,2);
delta_t=1e-4;T=10;t=0:delta_t:T;nt=size(t,2);
% system parameters
C=70.27;S0=1368;gama=0.67;theta=5.67e-8;
f=1/C*(1/4*(0.5+0.2*tanh((z-265)/10))*S0-gama*theta*z.^4);
d=0.0001;epsilon=5;alpha=1.9;

% initial conditions
Mf=max(abs(f));
PA0=sqrt(40/pi)*exp(-40*(z-228).^2);
PB0=sqrt(40/pi)*exp(-40*(z-279.7).^2);

% coefficients of the discrete equation
Calpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*sqrt(pi)*gamma(1-alpha/2));
Ch=d/2*(2/(zmax-zmin))^2-(2*epsilon/(zmax-zmin))^alpha*Calpha*zeta(alpha-1)*h^(2-alpha);
A=-2/(zmax-zmin);
B=-Calpha/alpha*(2*epsilon/(zmax-zmin)).^alpha*(1./(1+x).^alpha+1./(1-x).^alpha);
C=Ch/h^2;
D=Calpha*(2*epsilon/(zmax-zmin))^alpha*h^(-alpha);

% relative coefficients of the sum of the Toeplitz matrix
a=1:2*I-2;a=a.^(-1-alpha);
To=toeplitz([0 a],[0 a]);
diagg=sum(To)+[a (2*I-1)^(-1-alpha)]+flip([a (2*I-1)^(-1-alpha)]);
zz=2*(2*I-1)*real(ifft([0 a 0 flip(a)]));

% third-order Runge-Kutta integral
PA=zeros(nt,nx);PB=zeros(nt,nx);
PA(1,:)=PA0;PB(nt,:)=PB0;
for i=1:nt-1
    PA1=PA(i,:)+delta_t*One_Dimensional_FPK(PA(i,:));
    PA2=3/4*PA(i,:)+1/4*PA1+1/4*delta_t*One_Dimensional_FPK(PA1);
    PA(i+1,:)=1/3*PA(i,:)+2/3*PA2+2/3*delta_t*One_Dimensional_FPK(PA2);
    PB1=PB(nt-i+1,:)+delta_t*One_Dimensional_conj_FPK(PB(nt-i+1,:));
    PB2=3/4*PB(nt-i+1,:)+1/4*PB1+1/4*delta_t*One_Dimensional_conj_FPK(PB1);
    PB(nt-i,:)=1/3*PB(nt-i+1,:)+2/3*PB2+2/3*delta_t*One_Dimensional_conj_FPK(PB2);
end

% drawing
k=zeros(1,nt);
PP=PA.*PB;
for i=1:nt
    p=PP(i,:);
    [M,J]=max(p);
    k(i)=z(J);
end
figure;
plot(t,k); 

