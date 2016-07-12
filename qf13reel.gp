\\ matsupplement


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
qfsolve13(A,B)=
local(n,tmp,P,c,AA,BB,PP,AAA,BBB,y,z,z2,zz,Q,QQ,r,v,w,R,M,tmp1,s,tmp2,E,V,Y,RR,S,m,mm,C);
n=#A;
tmp=qfsign(A);
if(min(tmp[1],tmp[2])<floor(n/2),
  a=dichoqfsign(A,B);
  A=a*A+B);
P=Redqf(A)[1];
A=P*A*P~;
B=P*B*P~;
c=sign(B[1,1]);
print(c);
AA=vecextract(A,Str(3".."n),Str(3".."n));
BB=vecextract(B,Str(3".."n),Str(3".."n));
PP=DoubleWitt(AA,BB);
AAA=PP*AA*PP~;
BBB=PP*BB*PP~;
y=PosNeg(AAA,BBB,c);
z=y*PP;
print("A ", z*AA*z~);
print("B ",z*BB*z~);
Q=BaseWitt(AA)[1];
tmp=Q*AA*Q~;
Q[2,]=Q[2,]/(tmp[1,2]);
QQ=Mordell2(Q*AA*Q~);
Q=QQ*Q;
AAA=Q*AA*Q~;
BBB=Q*BB*Q~;
z2=bestappr(z,10^50);
zz=z2*(Q^-1); 
s=sign(eval(zz*BBB*zz~));
print(s);
while(s==c ,r=r+10; z2=bestapprox(z,10^100);zz=z2*(Q^-1);; s=sign(eval(zz*BBB*zz~));print(s);); 
\\print(z*AA*z~==zz*AAA*zz~);
\\print(z*(Q^-1*AAA*(Q^-1)~)*z~==(z*(Q^-1))*AAA*((Q^-1)~*z~));
\\print(z*AA*z~==z*(Q^-1*AAA*(Q^-1)~)*z~);
\\r=50;
v=Final(AAA,BBB,zz);
\\print("PosNeg :",v*AAA*v~," ",v*BBB*v~);
\\s=sign(eval(real(v*BBB*v~)));
\\while(s==c ,r=r+10; v=Final(AAA,BBB,zz,r); s=sign(eval(real(v*BBB*v~))); print("PosNeg :",v*AAA*v~," ",v*BBB*v~," ","r= ",r)); 
\\print("PosNeg :",v*AAA*v~," ",v*BBB*v~," ","r= ",r);
w=v*Q;
w=w/content(w);
\\print("w ",w);
R=matid(n-2);
\\R[1,]=w;
\\R[2]=vecteur(w);
M=Mat(w);
R=(mathnf(M,4)[2])^-1;
tmp1=R[1,];R[1,]=R[n-2,];R[n-2,]=tmp1;
E=R*AA*R~;
M=matrix(2,n-2);
M[1,1]=1;
\\V=M[1,]*E;
Y=vecteur(M[1,]*E);
M[2,]=Y;
RR=(mathnf(M,4)[2])^-1;
tmp2=RR[1,];RR[1,]=RR[n-3,];RR[n-3,]=tmp2;
tmp2=RR[2,];RR[2,]=RR[n-2,];RR[n-2,]=tmp2;
R=RR*R;
R[2,]=R[2,]/(R*AA*R~)[1,2];
AA=R*AA*R~;
RR=Mordell2(AA);
R=RR*R;
AA=RR*AA*RR~;
RR=Mordell3(AA);
R=RR*R;
S=matid(n);
for(i=3,n,for(j=3,n, S[i,j]=R[i-2,j-2]));
Q=S*P;


AA=S*A*S~;
BB=S*B*S~;

\\AAA=vecextract(AA,Str(5".."n),Str(5".."n));
\\BBB=vecextract(BB,Str(5".."n),Str(5".."n));

\\[r,s]=qfsign(AA);
\\u=min(r-2,s-2);
i=5;
\\t=-1;
P=matid(n);
while((n-i+1)>4,        
          F=vecextract(AA,Str(i".."n),Str(i".."n));
          tmp2=Redqf(F);
          if(type(tmp2)=="t_INT",break);
          tmp2=tmp2[1];
          C=matid(n);
          for(j=i,n, for(k=i,n, C[j,k]=tmp2[j-i+1,k-i+1]) );
          P=C*P; 
          AA=C*AA*C~;
          i=i+2;t=t+1
    );
R=P;
\\AA=R*A*R~;
\\BB=R*B*R~;
BB=R*BB*R~;

for(i=1,n,
          if(AA[i,i]!=0,m=i;break)
   );
m=m-1;
mm=m/2;

C=matrix(mm,mm);
for(i=1,mm-1,
          for(j=1,mm-1,C[i,j]=BB[2*i-1,2*j-1])
   );
for(i=1,mm-1,C[mm,i]=BB[m-1,2*i-1]);
for(i=1,mm-1,C[i,mm]=BB[2*i-1,m-1]);
C[mm,mm]=BB[m-1,m-1];
 
\\ X=qfsolve(C);     \\mettre un compteur temps!
ZZ=0;
while(ZZ==0 || type(X=qfsolve(ZZ*C*ZZ~))=="t_INT", ZZ=issmooth(C); print1("."));
 if(type(X=qfsolve(ZZ*C*ZZ~))!="t_INT",
    if(type(X)=="t_MAT", X=X[,1]);
   X=X~;
     X=X*ZZ);   \\,

\\print(X~*C*X);

Y=vector(n);
for(i=1,mm,Y[2*i-1]=X[i]);
return(Y*(R*Q));
}
















