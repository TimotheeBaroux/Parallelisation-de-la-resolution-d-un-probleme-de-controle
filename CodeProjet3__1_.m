function CodeProjet3
TempsParallele = 0;
tic;
NT=10; #nombre d'itérations de la résolution parallèle
N=50; #nombre de subdivisions de [0,T] pour la parallélisation
M=199; #discrétisation en temps
T=10;
Alpha=0;
A=[1 0 ; 0 2];
B=[0 1 ;1 0];
I=[1 0 ; 0 1];
yo=[1;0];
yc=[0;1];
c=rand(1,N*M+1);

for e=1:NT 
[y01,p01]=ValInt(T,yo,yc,M,N,Alpha,c,B,A);    
ck=[]; 
TempPreliminaire = toc;
TempsMachine = 0;
for l=1:N
  tic;
  ti=(l-1)*T/N;
  tj=l*T/N;
  Li=(l-1)*M+1;
  Lj=l*M+1;
  ci=c(Li:Lj-1);
  ya=(1-(l-1)/N)*y01(:,l)+p01(:,l)*(l-1)/N; # =Lambda(L-1)
  yb=(1-l/N)*y01(:,l)+ p01(:,l)*l/N;
  [ybest,pbest,cbest]=Csolve(ti,tj,ya,yb,M,Alpha,B,A,ci);
  ck=[ck cbest];
  tps=toc;
  TempsMachine = max(TempsMachine, tps);
end
ybest(:,M+1)
ck=[ck rand(1,1)];  
c=ck;
end
TempsParallele = TempPreliminaire + TempsMachine + TempsParallele
#tic;
#[yL,pL] = Csolve(0,T,yo,yc,N*M,Alpha,B,A,c);
#Tempstotal=toc;   #temps sans parallélisation
#Tempsfullefficiency = Tempstotal/N
endfunction

function [yL,pL,cbest] = Csolve(ti,tj,yo,yc,M,Alpha,B,A,c) 
  #recherche du meilleur c sur [TL-1,TL]
Rho = 0.1;
maxiter =10;
[J,yL,pL,grad]=fonctionelle(ti,tj,yo,yc,M,Alpha,c,B,A);

for k=1:maxiter
  cbest=c-Rho*grad;
  [J,yL,pL,grad]=fonctionelle(ti,tj,yo,yc,M,Alpha,c,B,A);
end
endfunction

function y= yget(M,yo,A,B,dt,c)
  y=zeros(2,M+1);
 y(:,1)=yo;
for n = 1:M
  y(:,n+1)=expm(-I*A*dt)*expm(-I*B*c(n)*dt)*y(:,n);
end
endfunction

function p= pget(M,yc,A,B,dt,c)
  p=zeros(2,M+1);
p(:,M+1)= yc; 
for n=M:-1:1
  p(:,n)=expm(I*B*c(n)*dt)*expm(I*A*dt)*p(:,n+1);
end
endfunction

function [J,y,p,grad]=fonctionelle(ti,tj,yo,yc,M,Alpha,c,B,A) 
  #calcul du gradient
t=linspace(ti,tj,M+1);
dt=t(2)-t(1);
y=yget(M,yo,A,B,dt,c);
p=pget(M,yc,A,B,dt,c);

J=-real(y(:,M+1)'*yc)+(Alpha/2)*dt*(c*c');
grad=zeros(1,M);
for n=1:M 
  grad(n)=Alpha*c(n) + real(p(:,n+1)'*(expm(-I*A*dt)*(I*B)*expm(-I*B*c(n)*dt)*y(:,n)));
end
endfunction

function [y,p]= ValInt(T,yo,yc,M,N,Alpha,c,B,A) 
  #génération des Lambda
t=linspace(0,T,N+1);
dt=t(2)-t(1);
y=zeros(2,N+1);
y(:,1)=yo;
for l = 1:N
  y(:,l+1)=expm(-I*A*dt)*expm(-I*B*c((l-1)*M+1)*dt)*y(:,l);
  for n = 1:M-1
    y(:,l+1)=expm(-I*A*dt)*expm(-I*B*c((l-1)*M+1+n)*dt)*y(:,l+1);
  end
end

p=zeros(2,N+1);
p(:,N+1)= yc; 
for l=N:-1:1
  p(:,l)=expm(I*B*c(l*M+1)*dt)*expm(I*A*dt)*p(:,l+1);
  for n=1:M-1
    p(:,l)=expm(I*B*c(l*M+1-n)*dt)*expm(I*A*dt)*p(:,l);
  end
end
endfunction