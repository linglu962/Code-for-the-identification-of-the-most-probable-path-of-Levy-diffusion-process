function AA=One_Dimensional_conj_FPK(P)
global f h I diagg zz A B C D;
% second-order derivatives
Ps1=[0 P];Ps1(end)=[];
Pp1=[P 0];Pp1(1)=[];

% calculating sum of the remaining terms using Toeplitz structure
PP=[P zeros(1,2*I-1)];
PP=real(fft(zz.*ifft(PP)));
PP=PP(1:(2*I-1));

% the whole equation
AA=-A*f.*(Pp1-Ps1)/2/h+B.*P+C*(Ps1-2*P+Pp1)+D*(PP-diagg.*P);

end