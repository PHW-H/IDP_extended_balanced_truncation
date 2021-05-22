function [H,R,J,B] = Modeltype42(Rl,Rc,Ll,Cc)
R=[Rl, Rc];

n=size(Cc);
n=n(1,1);

B=zeros(2*n,1);
B(n+1,1)=1/R(1,2);

l=size(Ll);
l=l(1,1);
rl=size(Rl);
rl=rl(1,1);
rc=size(Rc);
rc=rc(1,1); 

p=min(Rc);
if n~=c && n~=rl && n~=rc 
    'dimentions do not match'
    return
end
 %% creating matrix F
 A11(1,1)=-R(1,1);
 A22(n,n)=-1/R(n,2);
 A22(n-1,n)=1/R(n,2);
 A22(n,n-1)=1/R(n,2);
 H(1,1)=1/Ll(1);
 H(1+n,1+n)=1/Cc(1);
 i=2;
 j=n-1;
 k=j;
 while i<=n
    A11(i,i)=-R(i,1);
    A22(j,j)=-1/R(j+1,2)-1/R(j,2);
    while k>1
    A22(j-1,j)=1/R(j,2);
    A22(j,j-1)=1/R(j,2);
    k=k-1;
    end
    H(i,i)=1/Ll(i);
    H(i+n,i+n)=1/Cc(i);
    j=n-i;
    k=j;
    i=i+1;
 end
 
 a=ones(1,n);
 A12=diag(a);
 A21=-1*A12.';
 
 O=zeros(n,n);
 R=[A11,O;O,A22];
 J=[O,A12;A21,O];
end

