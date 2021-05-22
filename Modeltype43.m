function [H,R,J,B] = Modeltype43(Rl,Rc,Ll,Cc)
R=[Rl, Rc];

n=size(Ll);
n=n(1,1);

B=zeros(2*n,1);
B(n+1)=1;

c=size(Cc);
c=c(1,1);
rl=size(Rl);
rl=rl(1,1);
rc=size(Rc);
rc=rc(1,1); 

if n~=c && n~=rl && n~=rc 
    'dimentions do not match'
    return
end
 %% creating matrix F
 i=1;
 while i<=n
    A11(i,i)=-R(i,1);
    if R(i,2)==0
        A22(i,i)=0;
    else
        A22(i,i)=-1/R(i,2);
    end
    H(i,i)=1/Ll(i);
    H(i+n,i+n)=1/Cc(i);
    i=i+1;
 end
 
 a=ones(1,n);
 b=ones(1,n-1);
 A121=diag(-a);
 A122=diag(b,-1);
 A12=A121+A122;
 O=zeros(n,n);

 A122=diag(b,-1);
 A21=-1*A12.';
 R=[A11,O;O,A22];
 J=[O,A12;A21,O];
end
 

