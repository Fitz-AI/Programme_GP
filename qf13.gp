/* FIXME utiliser Qfsolve de Castel pourles changements de variables 
   quand les coefficients sont trop grands */
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





























