
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ---- Model ---- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%input model MANUAL
% Rc=[]; %values for risistance in Rc (in Ohm)
% Rl=[]; %values for risistance in Rl (in Ohm)
% Cc=[]; %values for capasitance in C (in Farad)
% Ll=[]; %values for inductance in L (in Henry)
% ModT=[]; %Model type
% M=[]; %reduction (in percentages)
% n=size(L);
% n=n(1,1);

%generat random model
% n=20; %dimentions of the system (number of conductors and capasitors)
% setMT=1; %determine model type or if=0 generates random model type
% setRe=0.25; %determine reduction in % or if=0 generates random reduction in %
% saveM=0; % if=0 saves model in 'Model_auto_save'
% [Rl,Rc,Cc,Ll,ModT,M] = Random_model_generator(n,setMT,setRe,saveM);

%load existing model
load ('Model_auto_save');

%%
if ModT==1
    [H,R,J,B] = Modeltype41(Rl,Rc,Ll,Cc);
elseif ModT==2
    [H,R,J,B] = Modeltype42(Rl,Rc,Ll,Cc);
elseif ModT==3
    [H,R,J,B] = Modeltype43(Rl,Rc,Ll,Cc);
else
    'error model type does not exist'
    return
end
F  =J+R;
A  =F*H;
Hi =inv(H);

beta=1; 
M   =2*M*n;

C  = B'*H;

if ModT==2
    do=0.10;
    Q = do*H;
    Xo = -Q*A-A'*Q-C'*C;
    EIGXo=eig(Xo);
    i=max(EIGXo);
    while i<0
        do=do*1.5;
        Q = do*H;
        C= B'*H;
        Xo = -Q*A-A'*Q-C'*C;
        EIGXo=eig(Xo);
        i=max(EIGXo);
    end
else
    syms 'do';
    Q  = do*H;
    Xo = -Q*A-A'*Q-C'*C; 

    EIGXo = eig(Xo);
    i = 1;
    while i<=2*n % find definition delta_o
        EIGXo(i)=solve(EIGXo(i)==0,do);
        i=i+1;
    end
    do = double(max(EIGXo))+0.0001;
end
Q  = do*H;
Qi=inv(Q);
Xo = -Q*A-A'*Q-C'*C; 

dc = do;

Pi = dc*Hi;
Xc = -A*Pi-Pi*A'-(B*B');

if Xc>=0 % find definition delta_c if needed
    syms 'dc'
    Pi=dc*Hi;
    Xc=-A*Pi-Pi*A'-(B*B');
    EIGXc=eig(Xc);
    i=1;
    while i<=2*n
        EIGXc(i)=solve(EIGXc(i)==0,dc);
        i=i+1;
    end
    dc=double(max(EIGXc))+0.0001;
    Pi=dc*Hi;
end

epsc = 1;
epso = 1;

i=1;
while i<=n %define zeta
    zeta_c(i,i)=Pi(i,i)*(1.1^i);
    zeta_o(i,i)=Q(i,i)*(1.1^i);
    zeta_c(n+i,n+i)=Pi(n+i,n+i)*(1.1^i);
    zeta_o(n+i,n+i)=Q(n+i,n+i)*(1.1^i);
    i=i+1;
end

GAMc=-epsc*zeta_c;
GAMo=zeta_o;

Thc=(-GAMc+A*Pi+B*B')*inv(Xc)*(-GAMc+Pi*A'+B*B');
condc=2*(beta*Pi+GAMc)-Thc;
con1=min(eig(condc)); % checking (all the eigenvalues must be positive)
i=1;
while con1<=0 %find smalles possible beta
    beta=beta*10;
    condc=2*(beta*Pi+GAMc)-Thc;
    con1=min(eig(condc));
    i=i+1;
    if i>11; %set maximum value for beta (if beta to large -> rounding errors)
        'error_beta'
        return
    end
end

epsc1=epsc;

i=1;
while con1>=0 %obtaining max value for epsilon_c
    epsc=epsc1;
    epsc1=(5*i);
    GAMc=-epsc1*zeta_c;
    Thc=(-GAMc+A*Pi+B*B')*inv(Xc)*(-GAMc+Pi*A'+B*B');
    condc=2*(beta*Pi+GAMc)-Thc;
    con1=min(eig(condc));
    i=i+1;
    if epsc1>=beta %epsc has to be smaller then beta
        'error_epsc'
        return
    end
end

GAMc=-epsc*zeta_c;

alpha=beta;

Tho=(GAMo-Q*A)*inv(Xo)*(GAMo-A'*Q); 

condo=2*(alpha*Q+GAMo)-Tho; 
con2=min(eig(condo)); %checking (all the eigenvalues must be positive)

epso1=epso;

i=1;
while con2 >= 0 
    epso=epso1;
    epso1=(5*i);
    GAMo=epso1*zeta_o;
    Tho=(GAMo-Q*A)*inv(Xo)*(GAMo-A'*Q);
    condo=2*(alpha*Q+GAMo)-Tho;
    con2=min(eig(condo));
    i=i+1;
    if epso>=beta %epsc has to be smaller then beta
        'error_epso'
        return
    end
end
GAMo=epso*zeta_o;

%%
%Define Ti

Ti=beta*Pi+GAMc;

min(eig(Ti)); % checking (all the eigenvalues must be positive)

Thc=(-GAMc+A*Pi+B*B')*inv(Xc)*(-GAMc+Pi*A'+B*B');
condc=2*(Ti)-Thc;
con1=min(eig(condc)); % checking (all the eigenvalues must be positive)
if con1<=0
    con1
    'error1'
    return
end
    

%%

% Define S

S=inv(alpha*Qi+Qi*GAMo*Qi);

Tho=(GAMo-Q*A)*inv(Xo)*(GAMo-A'*Q);

condo=2*(alpha*Q+GAMo)-Tho; 
con2=min(eig(condo)); %checking (all the eigenvalues must be positive)
if con2<=0
    con2
    'error2'
    return
end
%%%%%%%%%%%%%%%%%
%% Transformation Extended

%splitting Ti and S in C and L related parts
TiL=Ti(1:n,1:n);
TiC=Ti(n+1:2*n,n+1:2*n);
SL=S(1:n,1:n);
SC=S(n+1:2*n,n+1:2*n);

PhTiC=chol(TiC);
[UTSC,S2TSC]=svd(PhTiC*SC*PhTiC');
STSC=sqrt(S2TSC);
WEC=PhTiC'*UTSC*sqrt(inv(STSC));
WECi=inv(WEC);

PhTiL=chol(TiL);
[UTSL,S2TSL]=svd(PhTiL*SL*PhTiL');
STSL=sqrt(S2TSL);
WEL=PhTiL'*UTSL*sqrt(inv(STSL));
WELi=inv(WEL);



WE=[WEL zeros(n,n); zeros(n,n) WEC];
WEi=inv(WE);

%% Transformation Generalized

PiL=Pi(1:n,1:n);
PiC=Pi(n+1:2*n,n+1:2*n);
QL=Q(1:n,1:n);
QC=Q(n+1:2*n,n+1:2*n);

PhPiC=chol(PiC);
[UQPC,S2QPC]=svd(PhPiC*QC*PhPiC');
SQPC=sqrt(S2QPC);
WGC=PhPiC'*UQPC*sqrt(inv(SQPC));
WGCi=inv(WGC);

PhPiL=chol(PiL);
[UQPL,S2QPL]=svd(PhPiL*QL*PhPiL');
SQPL=sqrt(S2QPL);
WGL=PhPiL'*UQPL*sqrt(inv(SQPL));
WGLi=inv(WGL);


WG=[WGL zeros(n,n); zeros(n,n) WGC];
WGi=inv(WG);


%%

e=@(k,n) [zeros(k-1,1);1;zeros(n-k,1)];
% M is the number of state you want to truncate
M=M;
K=(n*2)-M;
aux1=[e(1,2*n)];
aux2=[e(n+1,2*n)];
i=2;
while i<=0.5*K
    aux3=[e(i,2*n)];
    aux4=[e(n+i,2*n)];
    aux1=[aux1,aux3];
    aux2=[aux2,aux4];
    i=i+1;
end
aux=[aux1,aux2];


% Reduced via extended

Ah=WE\A*WE;
Hh=WE'*H*WE;
Bh=WE\B;
C=B'*H;
Ch=C*WE;
Ar=aux'*Ah*aux;
Br=aux'*Bh;
Cr=Ch*aux;
Hr=aux'*Hh*aux;

% Reduced via generalized

Ahg=WG\A*WG;
Hhg=WG'*H*WG;
Bhg=WG\B;
Chg=C*WG;
Arg=aux'*Ahg*aux;
Brg=aux'*Bhg;
Crg=Chg*aux;
Hrg=aux'*Hhg*aux;

%%

lCn = eig(STSC)/max(eig(STSC));
lLn = eig(STSL)/max(eig(STSL));

lCni = flip(lCn);
lLni = flip(lLn);

%plot eigenvalues
% figure 
% plot(lCni,'bO','LineWidth',2)
% grid on
% title('Eigenvalues of $\Lambda_{ST_{1}}$','Interpreter','latex')
% xticks([0:n])
% figure
% plot(lLni,'rO','LineWidth',2)
% grid on
% title('Eigenvalues of $\Lambda_{ST_{2}}$','Interpreter','latex')
% xticks([0:n])



%%

% Error system

Ae = [Ah zeros(2*n,2*n-M); zeros(2*n-M,2*n) Ar];
Be = [Bh; Br];
Ce = [Ch -Cr];

Aeg = [Ahg zeros(2*n,2*n-M); zeros(2*n-M,2*n) Arg];
Beg = [Bhg; Brg];
Ceg = [Chg -Crg];

%%

% H inf normst

% extended

fsys = ss(A,B,C,0);
bsys = ss(Ah,Bh,Ch,0);
rsys = ss(Ar,Br,Cr,0);
esys = ss(Ae,Be,Ce,0);

[ninff,fpeakf] = hinfnorm(fsys);
[ninfb,fpeakb] = hinfnorm(bsys);
[ninfr,fpeakr] = hinfnorm(rsys);
[ninfe,fpeake] = hinfnorm(esys);

bgsys = ss(Ahg,Bhg,Chg,0);
rgsys = ss(Arg,Brg,Crg,0);
egsys = ss(Aeg,Beg,Ceg,0);

[ninfbg,fpeakbg] = hinfnorm(bgsys);
[ninfrg,fpeakrg] = hinfnorm(rgsys);
[ninfeg,fpeakeg] = hinfnorm(egsys);

[ninfe;ninfeg]