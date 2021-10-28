function AA=One_Dimensional_FPK(P)
global f h Mf I diagg zz A B C D;
% third-order WENO scheme
fP_pos=1/2*(f.*P+Mf*P);fP_neg=1/2*(f.*P-Mf*P);
phi=fP_pos;
phip1=[phi 0];phip1(1)=[];
phis1=[0 phi];phis1(end)=[];
phis2=[0 phis1];phis2(end)=[];
r_neg=(10^(-6)+(phi-2*phis1+phis2).^2)./(10^(-6)+(phip1-2*phi+phis1).^2);omega_neg=(1+2*r_neg.^2).^(-1);
phix_neg=1/(2*h)*(phip1-phis1)-omega_neg/(2*h).*(3*phis1-phis2-3*phi+phip1);
fPx_pos=phix_neg;
 
phi=fP_neg;
phip1=[phi 0];phip1(1)=[];
phip2=[phip1 0];phip2(1)=[];
phis1=[0 phi];phis1(end)=[];
r_pos=(10^(-6)+(phip2-2*phip1+phi).^2)./(10^(-6)+(phip1-2*phi+phis1).^2);omega_pos=(1+2*r_pos.^2).^(-1);
phix_pos=1/(2*h)*(phip1-phis1)-omega_pos/(2*h).*(3*phi-phis1-3*phip1+phip2);
fPx_neg=phix_pos;

fPx=fPx_pos+fPx_neg;

% second-order derivatives
Ps1=[0 P];Ps1(end)=[];
Pp1=[P 0];Pp1(1)=[];

% calculating sum of the remaining terms using Toeplitz structure
PP=[P zeros(1,2*I-1)];
PP=real(fft(zz.*ifft(PP)));
PP=PP(1:(2*I-1));

% the whole equation
AA=A*fPx+B.*P+C*(Ps1-2*P+Pp1)+D*(PP-diagg.*P);
end