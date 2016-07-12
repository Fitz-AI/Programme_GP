
{
Amer(A,B,z)=
local(n,P,PP,BB,Q,QQ);
n=#A;
P=matid(n);
P[1,]=z;
P[2,]=P[2,]/(P*A*P~)[1,2];
PP=Mordell2(P*A*P~);
P=PP*P;
BB=P*B*P~;
Q=matid(n);
Q[3,]=Q[3,]/BB[1,3];
QQ=Mordell6(Q*BB*Q~);
Q=QQ*Q;
P=Q*P;
return(P);
}



{
permute(A)=
local(n,P,j,k);
n=#A;
P=matrix(n,n);
P[1,1]=1;
j=0;k=1;
for ( i=2 , n , 
               if ( i % 2 == 0, P[i,n-j]= 1 ; j=j+1 ,
                                                   P[i,1+k]=1 ; k=k+1);
              );
return(P~);
\\P=P~;
\\AA=P*A*P~;
\\r=m/2;
\\Q0=vecextract(AA,Str(r+1".."n-r),Str(r+1".."n-r));
\\PP=qflllgram_indef(Q0)[2]~;
\\QQ=matid(n);
\\for (i=r+1,n-r, for (j=r+1,n-r, QQ[i,j]=PP[i-r,j-r]));
\\return(QQ*P);
}

{
s(a,x,Q)=
local(r);
r=x-(2*(a*Q*x~)/(a*Q*a~))*a;
return(r);
}


{
SOn(A)=
local(II,S,a);
n=#A;
II=matid(n);
S=matrix(n,n);
a=vector(n);
while ( (a*A*a~)==0 , a=vector(n,i,random(4)-2));
for (i=1,n, S[i,]=s(a,II[i,],A));
return(S);
} 

{
isotrope(A)=
local(n,m,F,D,P,Q0,tmp,tmp2,r,U,E,J,C);
n=#A;
for(i=1,n,
          if(A[i,i]!=0,r=i;break)
   );
P=matid(n);
r=r-1;
m=2*r;
tmp=0;
J=matrix(r,r);
for (i=1,r, for(j=1,r, if (i+j==r+1, J[i,j]=1)));
Q0=vecextract(A,Str(r+1".."r+(n-m)),Str(r+1".."r+(n-m)));
\\while(tmp==0, F=matrix(r,r,i,j,random(10)-5); tmp=matdet(F));
\\D=J*((F~)^-1)*J~;
\\for (i = 1 , r, for ( j= 1,r, P[i,j]=F[i,j]));
\\for (i = m , n , for (j=m,n, P[i,j]=D[i-m+1,j-m+1]));
U=matrix(n-m,r,i,j,random(10000000)-5000000);
 for (i = r+1 , n-r , for (j=1 , r , P[i,j]=U[i-r,j]));
E=matrix(r,r,i,j, random(1000)-500); E=E-E~;
T=-J*U~*Q0;
for (i = n-r+1 , n , for (j=r+1 , n-r , P[i,j]=T[i-(n-r),j-r]));
C=(-1/2)*J*(U~*Q0*U)+J*E;
for (i = n-r+1 , n , for (j=1 , r , P[i,j]=C[i-(n-r),j]));
return(P~);
}

{
isotrope2(A)=
local(n,m,F,D,P,Q0,M0,tmp,tmp2,r,U,E,J,C);
n=#A;
for(i=1,n,
          if(A[i,i]!=0,r=i;break)
   );
P=matid(n);
r=r-1;
m=2*r;
tmp=0;
J=matrix(r,r);
for (i=1,r, for(j=1,r, if (i+j==r+1, J[i,j]=1)));
Q0=vecextract(A,Str(r+1".."r+(n-m)),Str(r+1".."r+(n-m)));
tmp=0;
\\while(tmp==0, F=matrix(r,r,i,j,random(10)-5); tmp=matdet(F));
\\D=J*((F~)^-1)*J~;
\\for (i = 1 , r, for ( j= 1,r, P[i,j]=F[i,j]));
\\for (i = m , n , for (j=m,n, P[i,j]=D[i-m+1,j-m+1]));
while(tmp==0, F=matrix(r,r,i,j,random(10)-5); tmp=matdet(F));
for (i = 1 , r , for (j=1 , r , P[i,j]=F[i,j]));
D=J~*((F~)^-1)*J;
for (i = n-r+1 , n , for (j=n-r+1 , n , P[i,j]=D[i-n+r,j-n+r]));
M0=matrix(n-m,n-m,i,j,random(20)-10);
for (i = r+1 , n-r , for (j=r+1 , n-r , P[i,j]=M0[i-r,j-r]));
U=matrix(n-m,r,i,j,random(20)-10);
 for (i = r+1 , n-r , for (j=1 , r , P[i,j]=U[i-r,j]));
E=matrix(r,r,i,j, random(20)-10); E=E-E~;
T=-J~*((F~)^-1)*U~*Q0*M0;
for (i = n-r+1 , n , for (j=r+1 , n-r , P[i,j]=T[i-(n-r),j-r]));
C=(-1/2)*J*((F~)^-1)*(U~*Q0*U)+J*((F~)^-1)*E;
for (i = n-r+1 , n , for (j=1 , r , P[i,j]=C[i-(n-r),j]));
return(P~);
}

{
Qfn(A,B)=
my(n,tmp,a,P,m,d,cpt,Q,C,ZZ,X,Y,cpt2);
\\signa(A,B);
n=#A;
tmp=qfsign(A);
if(min(tmp[1],tmp[2])<floor(n/2),
  a=dichoqfsign(A,B);
  A=a*A+B);
P=test1(A);
A1=P*A*P~;
B1=P*B*P~;
for(i=1,n,
          if(A1[i,i]!=0,m=i;break)
   );
m=m-1;
P1=permute(A1);
P=P1*P;
A2=P1*A1*P1~;
B2=P1*B1*P1~;
mm=m/2;
print(mm);
\\print("m=",m);
cpt=0;
Q=P; A0=A2; B0=B2;
while( d==0,
if (cpt>0,
          \\if (cpt > 2*n , cpt=1; P=Q ; A2=A0; B2=B0; cpt2=cpt2+1;printf("."));
          \\S=SOn(A2); P=P*S; A2=S*A2*S~; B2=S*B2*S~;);
          \\A2=A0; B2=B0; P=Q; 
          S=isotrope2(A2); P=S*P; A2=S*A2*S~; B2=S*B2*S~;);
                                                                                                                                                                                                                                                                                                                                                      C=vecextract(B2,Str(1".."mm),Str(1".."mm));
if (mm>5, print("gaffe"); mm=5; C=vecextract(C,"1..5","1..5"));

if( #C==5 , C=C/content(C);
          if( type(X=myqfsolve(C)) != "t_INT", d=1 ; print(X), print("p= ", X); cpt=cpt+1 ;) 
   );

\\if( mm > 4 , 
\\          if( type(X=qfsolve(C)) != "t_INT", d=1 ; print(X), print1(" p= ", X); cpt=cpt+1 ;) 
\\   );
\\ faire un test pour savoir s'il y a des solutions pour n=5,4,3,2.
if (mm < 5 , R=factor(matdet(C),0)[,1];

\\print(factor(matdet(C)*denominator(matdet(C)),0));
          \\if( ispseudoprime(abs(matdet(C)*denominator(matdet(C))))==1, d=1; type(X=qfsolve(C)) != "t_INT", cpt=cpt+1;)
          if(matdet(C)==1  && type(X=qfsolve(C)) != "t_INT", d=1 , cpt=cpt+1;);
          if(ispseudoprime(R[#R])==1  && type(X=qfsolve(C)) != "t_INT", d=1 , print1(".");cpt=cpt+1;)
             );
  ); 
print("compteur ", cpt );
if(type(X)=="t_MAT", X=X[,1]);
X=X~;
Y=vector(n);
for(i=1,mm,Y[i]=X[i]);
return(Y*P);
}
