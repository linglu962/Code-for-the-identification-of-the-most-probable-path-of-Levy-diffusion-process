% step length of the integration
T=2;delta_t=1e-2;
t=0:delta_t:T;
% system parameters
alpha=0.5;epsilon=1;d=1;
% initial condition
x0=-1;
% integrating
for j=1:10000
    x=zeros(1,size(t,2));x(1)=x0;
    for i=1:size(t,2)-1
        f=x(i)-x(i).^3;
        x(i+1)=x(i)+f*delta_t+sqrt(delta_t)*sqrt(d)*randn(1,1)+epsilon*delta_t^(1/alpha)*stblrnd(alpha,0,1,0,1,1);
    end
    if x(end)<1+1e-3&&x(end)>1-1e-3
        break
    end
end
plot(t,x);