function A=SSR(X0,Y0,k)
A=(X0'*X0)\(X0'*Y0);
X=X0;Y=Y0;c=A;p=1:size(A,1);
for i=1:k
    [MM,I]=min(abs(c));
    X(:,I)=[];
    c=(X'*X)\(X'*Y);
    p(abs(A)==MM)=0;
    A(abs(A)==MM)=0;
    A(find(p))=c; %#ok<FNDSB>
end