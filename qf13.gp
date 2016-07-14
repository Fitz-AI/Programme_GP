\r DoubleWitt.gp
\r qf2diag.gp
\r equilibre.gp
\r Qfmodife.gp
{
Mordell2(A)=
local(n,P);
n=#A;
P=matid(n);
for(i=3,n, P[i,2]=-A[i,1]);
return(P);
}

{
Mordell3(A)=
local(n,P,S);
n=#A;
P=matid(n);
for(i=2,n, P[i,1]=-A[i,2]);
A2=P*A*P~;
S=matid(n);
S[2,]=2*S[2,]-A2[2,2]*S[1,];
S[1,1]=1/2;
P=S*P;
return(P);
}

{
BaseWitt2(A)=
my(n,P,u,Z,X,Y,M,tmp1,p,tmp2,B);
n=#A;
P=matid(n);
A=A/content(A);
X=qfsolve2(A);\\("borne de Cassels= \n"round(vecmax(abs(A)^((n-1)/2))));
if(type(X)=="t_MAT", X=X[,1]~);
if(type(X)=="t_INT",return(1));
if(type(X)=="t_COL", X=(X/content(X))~);
\\M=Mat(X);
\\P=(mathnf(M,4)[2])^-1;
\\tmp1=P[1,];P[1,]=P[n,];P[n,]=tmp1;
Z=vector(n);
while( X*A*Z~==0, 
                 u=vector(n,i,random(50)-20); 
                 while (X*A*u~==0, u=vector(n,i,random(50)-20)); 
                 Z=u-((u*A*u~)/(2*(X*A*u~)))*X);
P=matsupplement(Mat([X~,Z~]))~;
\\P=matsupplement(Mat(X)~)~;
\\A2=P*A*P~;
\\M=matrix(2,n);
\\M[1,1]=1;
\\V=M[1,]*A2;
\\Y=vecteur(M[1,]*A2);
\\M[2,]=Y;
\\PP=(mathnf(M,4)[2])^-1;
\\tmp2=PP[1,];PP[1,]=PP[n-1,];PP[n-1,]=tmp2;
\\tmp2=PP[2,];PP[2,]=PP[n,];PP[n,]=tmp2;
\\P=PP*P;
return([P,P*A*P~]);
}

{
Redqf2(Q)=
local(P,QQ,PP,tmp,g);
tmp=BaseWitt2(Q);
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
Hyperbolique(A)=
my(r,u,s,n,tmp,B,P,tmp2,C);
n=#A;
[r,s]=qfsign(A);
u=min(r,s);
i=1;t=-1;
P=matid(n);
while(t<u,        
          B=vecextract(A,Str(i".."n),Str(i".."n));
          tmp2=Redqf2(B);
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
qf13(A,B)=
my(n,tmp,a,y,R,A1,B1,z1,z2,P,P2,A2,B2,z3,r,z,s,u,v,P3,P4,A3,A4,S2,A5,w,x);
n=#A;
tmp=qfsign(A);
if(min(tmp[1],tmp[2])<floor(n/2),
  a=dichoqfsign(A,B);
  A=a*A+B);
y=qfsolve(A);
if(type(y)=="t_MAT", y=y[,1]~);
if(type(y)=="t_COL", y=(y/content(y))~);
if ( y*B*y~<0, B=-B);
R=DoubleWitt(A,B);
A1=R*A*R~;
B1=R*B*R~;
z1=PosNeg(A1,B1,1);
z2=z1*R;
P=BaseWitt2(A)[1];
tmp=P*A*P~;
P[2,]=P[2,]/tmp[1,2];
P2=Mordell2(P*A*P~);
P=P2*P;
A2=P*A*P~;
B2=P*B*P~;
z3=z2*P^-1;
r=1;
z=Final(A2,B2,z3,r);
s=sign(eval(z*B2*z~));
while(s==1, r=2*r; z=Final(A2,B2,z3,r); s=sign(eval(z*B2*z~)));
z=z*P;
\\z=z/content(z);
while( y*A*z~==0, 
                 u=vector(n,i,random(50)-20); 
                 while (y*A*u~==0, u=vector(n,i,random(50)-20)); 
                 v=u-((u*A*u~)/(2*(y*A*u~)))*y; 
                 if (v*B*v~==0, return(v));
                 if (v*B*v~>0, y=v; z=v)
);
P3=matintersect(matker(Mat(y*A)),matker(Mat(z*A)))~;
P4=matid(n);
P4[1,]=y;
P4[2,]=z;
for(i=3,n, P4[i,]=P3[i-2,]);
A3=P4*A*P4~;
B3=P4*B*P4~;
C=A3[3..n,3..n];
D=B3[3..n,3..n];
C=C/content(C);
S=Hyperbolique(C);
P5=matid(n);
for(i=3,n,
          for( j=3,n,
                     P5[i,j]=S[i-2,j-2]));
A4=P5*A3*P5~;
B4=P5*B3*P5~;
P=P5*P4;
if ( B4[3,3]>0, 
                S2=matid(n);
                tmp=S2[1,];
                S2[1,]=S2[2,];
                S2[2,]=tmp;
                P=S2*P;
                A4=S2*A4*S2~;
                B4=S2*B4*S2~;
);
A5=matrix(5,5);
for(i=1,4,
          for(j=1,4,A5[i,j]=B4[2*i-1,2*j-1])
   );
for(i=1,4,A5[5,i]=B4[9,2*i-1]);
for(i=1,4,A5[i,5]=B4[2*i-1,9]);
A5[5,5]=B4[9,9];
w=myqfsolve(A5)~;
print(" Attention ", w*A5*w~);
x=vector(n);
for(i=1,5,x[2*i-1]=w[i]);
return(x*P);
}





























