clear;clc;
% mesh parameters
zmin=-2;zmax=2;M=1e7;
z=linspace(zmin,zmax,M);
delta_t=0.001;
% system parameter
f=z-z.^3;
d=0;epsilon=100;alpha=0.5;
% generating data
x=z+f*delta_t+delta_t^(1/2)*sqrt(d)*randn(1,M)+epsilon*delta_t^(1/alpha)*stblrnd(alpha,0,1,0,1,M);
R=abs(x-z);

% identification of stable parameter ¦Á,levy noise intensity ¦Å
N=2;e=1;m=5;
n=zeros(1,N+1);
for i=1:N+1
    n(i)=size(find(R>=e*m^(i-1)&R<e*m^i),2);
end
k=1:N;
alpha1=(k*log(m)).^(-1).*log(n(1)./n(2:end));
alpha_test=sum(alpha1)/N;

k=0:N;c=alpha_test*gamma((1+alpha_test)/2)/(2^(1-alpha_test)*sqrt(pi)*gamma(1-alpha_test/2));
epsilon1=(alpha_test*e^alpha_test*m.^(k*alpha_test).*n/2/c/delta_t/M/(1-m^(-alpha_test))).^(1/alpha_test);
epsilon_test=sum(epsilon1)/(N+1);

% identification of the drift term
p=find(R<e);
z1=z(p);x1=x(p);M1=size(z1,2);R1=x1-z1;
%binning
Q=50;% Number of bins
z1Q=zeros(1,Q);R1Q=zeros(1,Q);Weight=zeros(1,Q);RR1Q=zeros(1,Q);
h=(zmax-zmin)/Q;
for i=1:Q
    aa=find(z1>=zmin+(i-1)*h&z1<zmin+i*h);
    z1Q(i)=sum(z1(aa))/size(aa,2);
    R1Q(i)=sum(R1(aa))/size(aa,2);
    RR1Q(i)=sum(R1(aa).^2)/size(aa,2);
    Weight(i)=size(aa,2)/size(R1,2);
end

A=diag(Weight)*[ones(Q,1) z1Q' z1Q'.^2 z1Q'.^3 z1Q'.^4 z1Q'.^5 z1Q'.^6];
B=M1/M/delta_t*diag(Weight)*R1Q';

X0=A;Y0=B;
% initial solution
c0=(X0'*X0)\(X0'*Y0);
X=X0;Y=Y0;
% intermediate solution and record matrix
c=c0;
cc=zeros(size(c0,1),size(c0,1));
% pointer vector
p=1:size(c0,1);
p1=1:floor(1/5*Q);p2=floor(1/5*Q)+1:floor(2/5*Q);p3=floor(2/5*Q)+1:floor(3/5*Q);p4=floor(3/5*Q)+1:floor(4/5*Q);p5=floor(4/5*Q)+1:Q;
% vector of CV scores
delta=zeros(1,size(c0,1));
% sparsity enforcement
for i=1:size(A,2)
    [MM,I]=min(abs(c));
    X(:,I)=[];
    c=(X'*X)\(X'*Y);
    p(abs(c0)==MM)=0;
    c0(abs(c0)==MM)=0; 
    c0(find(p))=c;  %#ok<FNDSB>
    cc(:,i)=c0;
    % calculation of CV scores
    part1=Y0(p1)-X0(p1,:)*SSR(X0([p2 p3 p4 p5],:),Y0([p2 p3 p4 p5]),i);
    part2=Y0(p2)-X0(p2,:)*SSR(X0([p1 p3 p4 p5],:),Y0([p1 p3 p4 p5]),i);
    part3=Y0(p3)-X0(p3,:)*SSR(X0([p1 p2 p4 p5],:),Y0([p1 p2 p4 p5]),i);
    part4=Y0(p4)-X0(p4,:)*SSR(X0([p1 p2 p3 p5],:),Y0([p1 p2 p3 p5]),i);
    part5=Y0(p5)-X0(p5,:)*SSR(X0([p1 p2 p3 p4],:),Y0([p1 p2 p3 p4]),i);
    delta(end-i+1)=sqrt(1/5*(sum(part1.^2)+sum(part2.^2)+sum(part3.^2)+sum(part4.^2)+sum(part5.^2)));
end
ratio=delta(1:size(A,2)-1)./delta(2:size(A,2));
p=find(ratio<1.1,1);
drift_term_test=cc(:,end-p+1); %#ok<NASGU>

% identification of the diffusion term
c=alpha_test*gamma((1+alpha_test)/2)/(2^(1-alpha_test)*sqrt(pi)*gamma(1-alpha_test/2));
S=2*epsilon_test^alpha_test*e^(2-alpha_test)*c*(2-alpha_test)^(-1);
B1=diag(Weight)*(M1/M/delta_t*RR1Q'-S);

X0=A;Y0=B1;

c0=(X0'*X0)\(X0'*Y0);
X=X0;Y=Y0;

c=c0;
cc=zeros(size(c0,1),size(c0,1));

p=1:size(c0,1);
p1=1:floor(1/5*Q);p2=floor(1/5*Q)+1:floor(2/5*Q);p3=floor(2/5*Q)+1:floor(3/5*Q);p4=floor(3/5*Q)+1:floor(4/5*Q);p5=floor(4/5*Q)+1:Q;

delta=zeros(1,size(c0,1));

for i=1:size(A,2)
    [MM,I]=min(abs(c));
    X(:,I)=[];
    c=(X'*X)\(X'*Y);
    p(abs(c0)==MM)=0;
    c0(abs(c0)==MM)=0; 
    c0(find(p))=c;  %#ok<FNDSB>
    cc(:,i)=c0;
    
    part1=Y0(p1)-X0(p1,:)*SSR(X0([p2 p3 p4 p5],:),Y0([p2 p3 p4 p5]),i);
    part2=Y0(p2)-X0(p2,:)*SSR(X0([p1 p3 p4 p5],:),Y0([p1 p3 p4 p5]),i);
    part3=Y0(p3)-X0(p3,:)*SSR(X0([p1 p2 p4 p5],:),Y0([p1 p2 p4 p5]),i);
    part4=Y0(p4)-X0(p4,:)*SSR(X0([p1 p2 p3 p5],:),Y0([p1 p2 p3 p5]),i);
    part5=Y0(p5)-X0(p5,:)*SSR(X0([p1 p2 p3 p4],:),Y0([p1 p2 p3 p4]),i);
    delta(end-i+1)=sqrt(1/5*(sum(part1.^2)+sum(part2.^2)+sum(part3.^2)+sum(part4.^2)+sum(part5.^2)));
end
ratio=delta(1:size(A,2)-1)./delta(2:size(A,2));
p=find(ratio<1.1,1);
d_test=cc(:,end-p+1);