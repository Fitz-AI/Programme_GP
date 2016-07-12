{
issmooth(f)=
local(nn,P,Q,N,R);
nn=#f;
P=matrix(nn-1,nn,i,j,random(10)-5);
Q=P*f*P~;
N=matdet(Q); if(N==0, return(0));
\\print(log(abs(N))/log(10));
R=factor(N,0)[,1];
if(ispseudoprime(R[#R])==1,return(P),return(0));
}

{
vecteur(x)=
local(n,Y,Z,u,v,tmp,a,b);
n=#x;
Y=vector(n);
Z=qflll(Mat(x));
for(i=1,n,Y[i]=Z[i,1]);
return(Y);
}

{
Mordell2(Q)=
local(n,P);
n=#Q;
P=matid(n);
for(i=3,n, P[i,2]=-Q[i,1]);
return(P);
}

{
Mordell3(Q)=
local(n,P,S);
n=#Q;
P=matid(n);
for(i=2,n, P[i,1]=-Q[i,2]);
Q=P*Q*P~;
S=matid(n);
S[2,]=2*S[2,]-Q[2,2]*S[1,];
S[1,1]=1/2;
P=S*P;
return(P);
}

{
BaseWitt(A)=
local(n,P,ZZ,X,Y,M,tmp1,p,tmp2,B);
n=#A;
P=matid(n);
X=qfsolve(A);\\("borne de Cassels= \n"round(vecmax(abs(A)^((n-1)/2))));
\\print(X);
\\while(ZZ==0,ZZ=issmooth(C);print1(":"));
\\X=qfsolve(ZZ*A*ZZ~)~;
\\X=X*ZZ;
if(type(X)=="t_MAT", X=X[,1]~);
if(type(X)=="t_INT",return(1));
if(type(X)=="t_COL", X=(X/content(X))~);
M=Mat(X);
P=(mathnf(M,4)[2])^-1;
tmp1=P[1,];P[1,]=P[n,];P[n,]=tmp1;
B=P*A*P~;
M=matrix(2,n);
M[1,1]=1;
V=M[1,]*B;
Y=vecteur(M[1,]*B);
M[2,]=Y;
PP=(mathnf(M,4)[2])^-1;
tmp2=PP[1,];PP[1,]=PP[n-1,];PP[n-1,]=tmp2;
tmp2=PP[2,];PP[2,]=PP[n,];PP[n,]=tmp2;
P=PP*P;
return([P,P*A*P~]);
}

{
Redqf(Q)=
local(P,QQ,PP,tmp,g);
tmp=BaseWitt(Q);
if(type(tmp)=="t_INT",return(1));
P=tmp[1]; 
QQ=tmp[2]; 
g=QQ[1,2];
QQ[1,]=QQ[1,]/g;QQ[,1]=QQ[,1]/g;P[1,]/=g;
PP=Mordell2(QQ);
QQ=PP*QQ*PP~ ;
P=PP*P;
PP=Mordell3(QQ);
P=PP*P;
tmp=P*Q*P~;
\\if(tmp[2,1]!=1,P[2,]=P[2,]/tmp[2,1]);
return([P,P*Q*P~]);
}

{
test1(A)=
local(r,s,n,tmp,B,P,tmp2,C);
n=#A;
[r,s]=qfsign(A);
u=min(r,s);
i=1;t=-1;
P=matid(n);
while(t<u,        
          B=vecextract(A,Str(i".."n),Str(i".."n));
          tmp2=Redqf(B);
          if(type(tmp2)=="t_INT",break);
          tmp2=tmp2[1];
          C=matid(n);
          for(j=i,n, for(k=i,n, C[j,k]=tmp2[j-i+1,k-i+1]) );
          P=C*P; 
          A=C*A*C~;
          i=i+2;t=t+1
    );
return(P);
}






{
BorneCauchy(P)=
local(Pn,v);
Pn=poldegree(P);
v=vector(Pn);
for(i=1,Pn,v[i]=abs(polcoeff(P,i)));
return(1+vecmax(v));
}

{
dichoqfsign(A,B)=
local(n,P,b,m,q,c,d);
n=#A;
\\if(min(qfsign(B)[1],qfsign(B)[2])==floor(n/2),return(0));
\\ test de m=1
if(min(qfsign(A+B)[1],qfsign(A+B)[2])==floor(n/2),return(1));
\\ test de -1
if(min(qfsign(-A+B)[1],qfsign(-A+B)[2])==floor(n/2),return(-1));
if(min(qfsign(2*A+B)[1],qfsign(2*A+B)[2])==floor(n/2),return(2));
if(min(qfsign(-2*A+B)[1],qfsign(-2*A+B)[2])==floor(n/2),return(-2));
P=matdet(A*x+B);
b=BorneCauchy(P);
a=-b;

\\test en +oo
if((qfsign(A)[1])==ceil(n/2), if(qfsign(10*A+B)[1]==ceil(n/2), return(10));   
   while(polsturm(P,b/10,b+1)==0,b=b/10);
\\ prendre la valeur entière   
   if(qfsign(floor(b)*A+B)[1]==ceil(n/2), return(floor(b))); \\ prendre la valeur entière
   if(qfsign(ceil(b)*A+B)[1]==ceil(n/2), return(ceil(b))); 
   return(b));

\\test en -oo
if((qfsign(A)[2])==ceil(n/2), if(qfsign(-10*A+B)[1]==ceil(n/2), return(-10));
   while(polsturm(P,-b-1,-b/10)==0,b=b/10);
\\prendre valeur entière
   if(qfsign(floor(-b)*A+B)[1]==ceil(n/2), return(floor(-b))); 
   if(qfsign(ceil(-b)*A+B)[1]==ceil(n/2), return(ceil(-b)));  
   return(-b));

II=[a,b];
m=1;
\\ réduire l'intervalle 
if((ceil(n/2)-qfsign(B)[1])*(ceil(n/2)-qfsign(100*A+B)[1])<0, a=0; b=100; m=50);
if((ceil(n/2)-qfsign(B)[1])*(ceil(n/2)-qfsign(-100*A+B)[1])<0, a=-100; b=0;m=-50);
if((ceil(n/2)-qfsign(B)[1])*(ceil(n/2)-qfsign(2*A+B)[1])<0, a=0; b=2; m=1);
if((ceil(n/2)-qfsign(B)[1])*(ceil(n/2)-qfsign(-2*A+B)[1])<0, a=-2; b=0;m=-1);

\\dichotomie
q=qfsign(A+B)[1];
while(m==0 || q!=ceil(n/2),
     c=qfsign(a*A+B)[1];
     d=qfsign(m*A+B)[1];
     if((ceil(n/2)-c)*(ceil(n/2)-d)<0, b=m, a=m);
     m=(a+b)/2;
     q=qfsign(m*A+B)[1]
     );
if(qfsign(floor(m)*A+B)[1]==ceil(n/2) && m!=0, return(floor(m))); 
if(qfsign(ceil(m)*A+B)[1]==ceil(n/2) && m!=0 , return(ceil(m)));   
return(m);
}
   
{
finalqf13(A,B)=
local(A1,B1,n,tmp,tmp2,tmp3,tmp4,a,P,m,cpt,d,AA,BB,u,AAA,BBB,v,AAAA,BBBB,w,C,ZZ,X,Y);
\\signa(A,B);
n=#A;
tmp=qfsign(A);
if(min(tmp[1],tmp[2])<floor(n/2),
  a=dichoqfsign(A,B);
  if( a!=0 ,A=a*A+B));
P=test1(A);
A1=P*A*P~;
B1=P*B*P~;
for(i=1,n,
          if(A1[i,i]!=0,m=i;break)
   );
m=m-1;
mm=m/2;
\\print("m=",m);

AA=vecextract(A1,Str(1".."m),Str(1".."m));
BB=vecextract(B1,Str(1".."m),Str(1".."m));
tmp2=qfsign(BB);


if(min(tmp2[1],tmp2[2])!=floor(#BB/2),     \\si A dégénéré?????
   u=dichoqfsign(AA,BB);
   BB=u*AA+BB);




AAA=vecextract(AA,Str(1".."m-1),Str(1".."m-1));
BBB=vecextract(BB,Str(1".."m-1),Str(1".."m-1));





AAAA=matrix(m-2,m-2);
for(i=1,m-3,
          for(j=1,m-3,AAAA[i,j]=AAA[i,j])
   );
for(i=1,m-3,AAAA[m-2,i]=AAA[m-1,i]);
for(i=1,m-3,AAAA[i,m-2]=AAA[i,m-1]);
AAAA[m-2,m-2]=AAA[m-1,m-1];

BBBB=matrix(m-2,m-2);
for(i=1,m-3,
          for(j=1,m-3,BBBB[i,j]=BBB[i,j])
   );
for(i=1,m-3,BBBB[m-2,i]=BBB[m-1,i]);
for(i=1,m-3,BBBB[i,m-2]=BBB[i,m-1]);
BBBB[m-2,m-2]=BBB[m-1,m-1];



C=matrix(mm,mm);
for(i=1,mm-1,
          for(j=1,mm-1,C[i,j]=BBBB[2*i-1,2*j-1])
   );
for(i=1,mm-1,C[mm,i]=BBBB[m-2,2*i-1]);
for(i=1,mm-1,C[i,mm]=BBBB[2*i-1,m-2]);
C[mm,mm]=BBBB[m-2,m-2];

if ( C==matrix(mm,mm) , return("erreur"));
tmp=qfsign(C);
if(tmp[1]==0 || tmp[2]==0 , print("problème de signature");return(qfsolve13(A,B)));


\\X=qfsolve(C);     \\mettre un compteur temps!
ZZ=0;
while(ZZ==0 || type(X=qfsolve(ZZ*C*ZZ~))=="t_INT", ZZ=issmooth(C); print1("."));
 if(type(X=qfsolve(ZZ*C*ZZ~))!="t_INT",
    if(type(X)=="t_MAT", X=X[,1]);
   X=X~;
     X=X*ZZ);   \\,

\\print(X~*C*X);

Y=vector(n);
for(i=1,mm,Y[2*i-1]=X[i]);
return(Y*P);
}
