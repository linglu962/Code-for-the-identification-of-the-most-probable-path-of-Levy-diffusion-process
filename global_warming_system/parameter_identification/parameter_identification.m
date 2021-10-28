clear;clc;
% mesh parameters
zmin=220;zmax=320;M=1e8;
z=linspace(zmin,zmax,M);
delta_t=0.001;
% system parameter
C=70.27;S0=1368;gama=0.67;theta=5.67e-8;
f=1/C*(1/4*(0.5+0.2*tanh((z-265)/10))*S0-gama*theta*z.^4);
alpha=0.5;epsilon=1;d=1;

% generating data and necessary parameters in learning
x=z+f*delta_t+delta_t^(1/2)*sqrt(d)*randn(1,M)+epsilon*delta_t^(1/alpha)*stblrnd(alpha,0,1,0,1,M);
n_basis_function=6;
Q=1000;% Number of bins

% coordinate transformation
zzmin=-1;zzmax=5;
zz=linspace(zzmin,zzmax,M);
xx=(zzmax-zzmin)/(zmax-zmin)*x+(zmax*zzmin-zmin*zzmax)/(zmax-zmin);
R=abs(xx-zz);

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
epsilon_test_p=sum(epsilon1)/(N+1);
epsilon_test=(zmax-zmin)/(zzmax-zzmin)*epsilon_test_p;

% binning
p=find(R<e);
zz1=zz(p);xx1=xx(p);M1=size(zz1,2);R1=xx1-zz1;
z1Q=zeros(1,Q);R1Q=zeros(1,Q);Weight=zeros(1,Q);RR1Q=zeros(1,Q);
h=(zzmax-zzmin)/Q;
for i=1:Q
    aa=find(zz1>=zzmin+(i-1)*h&zz1<zzmin+i*h);
    z1Q(i)=sum(zz1(aa))/size(aa,2);
    R1Q(i)=sum(R1(aa))/size(aa,2);
    RR1Q(i)=sum(R1(aa).^2)/size(aa,2);
    Weight(i)=size(aa,2)/size(R1,2);
end

A=diag(Weight)*basis_diff(z1Q,n_basis_function);
B=M1/M/delta_t*diag(Weight)*R1Q';

% identification of the diffusion term
c=alpha_test*gamma((1+alpha_test)/2)/(2^(1-alpha_test)*sqrt(pi)*gamma(1-alpha_test/2));
S=2*epsilon_test_p^alpha_test*e^(2-alpha_test)*c*(2-alpha_test)^(-1);
B1=diag(Weight)*(M1/M/delta_t*RR1Q'-S);

X0=A;Y0=B1;
% initial solution
c0=(X0'*X0)\(X0'*Y0);
X=X0;Y=Y0;

% intermediate solution and record matrix
c=c0;
cc=zeros(n_basis_function,n_basis_function);
cc(:,1)=c0;

% pointer vector
p=1:n_basis_function;

% sparsity enforcement
for i=1:n_basis_function-1
    [MM,I]=min(abs(c));
    X(:,I)=[];
    c=(X'*X)\(X'*Y);
    p(abs(c0)==MM)=0;
    c0(abs(c0)==MM)=0; 
    c0(find(p))=c;  %#ok<FNDSB>
    cc(:,i+1)=c0;
end
% Cross Validation
% run times and numbers of folds
Nk=100;fold=5;
% vector of CV scores
delta=zeros(1,n_basis_function);
% calculation of CV scores
for i=1:Nk
    %  data grouping
    s = randperm(Q);
    ngroup=zeros(1,fold);
    ngroup(1) = ceil(rand*(Q-fold+1));
    for j=2:fold-1
        ngroup(j) = ceil(rand*(Q-sum(ngroup(1:j-1))-fold+j));
    end
    ngroup(end)=Q-sum(ngroup(1:fold-1));
    for j=1:n_basis_function
        part=zeros(1,fold);
        p1=s(1:ngroup(1));
        p2=s;p2(1:ngroup(1))=[];
        part(1)=sum((Y0(p1)-X0(p1,:)*SSR(X0(p2,:),Y0(p2),j-1)).^2);
        for m=2:fold
            p1=s(sum(ngroup(1:m-1))+1:sum(ngroup(1:m)));
            p2=s;p2(sum(ngroup(1:m-1))+1:sum(ngroup(1:m)))=[];
            part(m)=sum((Y0(p1)-X0(p1,:)*SSR(X0(p2,:),Y0(p2),j-1)).^2);
        end
        delta0=1/fold*sum(part);
        delta(end-j+1)=delta(end-j+1)+delta0;
    end
end
delta=sqrt(delta/Nk);
ratio=delta(1:n_basis_function-1)./delta(2:n_basis_function);
% drawing of CV scores and ratios
figure
plot(1:n_basis_function-1,delta(1:end-1));
hold off
figure
plot(2:n_basis_function,ratio(1:end));
hold off


d_test=cc_true(:,end)';