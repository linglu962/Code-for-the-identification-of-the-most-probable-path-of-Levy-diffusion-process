% step length of the integration
T=10;delta_t=1e-2;
t=0:delta_t:T;
% system parameters
C=70.27;S0=1368;gamma=0.67;theta=5.67e-8;
alpha=0.5;epsilon=1;d=0.0001;
% initial condition
x0=228;
% integrating
for j=1:10000
    x=zeros(1,size(t,2));x(1)=x0;
    for i=1:size(t,2)-1
        f=1/C*(1/4*(0.5+0.2*tanh((x(i)-265)/10))*S0-gamma*theta*x(i).^4);
        x(i+1)=x(i)+f*delta_t+sqrt(delta_t)*sqrt(d)*randn(1,1)+epsilon*delta_t^(1/alpha)*stblrnd(alpha,0,1,0,1,1);
    end
    if x(end)<279.7+1e-1&&x(end)>279.7-1e-1
        break
    end
end
plot(t,x);