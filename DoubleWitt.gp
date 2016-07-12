{
Mordell2(Q)=
local(n,P);
n=#Q;
P=matid(n);
for(i=3,n, P[i,2]=-Q[i,1]/Q[2,1]);
return(P);
}

{
Mordell6(Q)=
local(n,P);
n=#Q;
P=matid(n);
P[2,3]=-Q[2,1];
for(i=4,n, P[i,3]=-Q[i,1]);
return(P);
}


{
DoubleWitt(A,B)=
local(n,z,P,PP,BB,Q,QQ);
n=#A;
\\z=qfreel(A,B);
z= finalsol(A,B);
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
PosNeg(A,B,c)=;
local(n,AA,d,BB,y,tmp,e,z);
n=#A;
AA=vecextract(A,Str(2".."n),Str(2".."n));
BB=vecextract(B,Str(2".."n),Str(2".."n));
d=sign(eval(real(AA[2,2])));
y=vector(n-1,i,1);
y[2]=(c*d);
tmp=(y*BB*y~)-y[2]*(y*AA*y~);
e=sign(eval(real(tmp)));
while ( e== 1 , y[2]=2*y[2] ; tmp= (y*BB*y~)-y[2]*(y*AA*y~);e=sign(eval(real(tmp))); print("e= "e);print(y); );
z=vector(n);
z[1]=-(y*AA*y~)/(2*y[1]);
for(i=2,n,z[i]=y[i-1]);
return(z);
}

\\{
\\Approx(v,r)=
\\local(n,M,T,z);
\\n=#v;
\\M=matid(n+1);
\\for(i=1,n, M[n+1,i]=v[i]);
\\M[n+1,n+1]=1/r;
\\T=qflll(M);
\\z=vector(n);
\\for(i=1,n,for(j=2,n+1,if (abs(abs(v[i])/(abs(T[n+1,j])/r))<(10/r),z[i]=T[n+1,j])));
\\for(i=1,n,if(sign(z[i])!= sign(v[i]),z[i]=-1*z[i]));
\\ return(z/r);
\\return(z);
\\}



{
Final2(A,B,z,r)=
local(n,AA,zz,y,yy);
n=#A;
AA=vecextract(A,Str(2".."n),Str(2".."n));
zz=vecextract(z,Str(2".."n));
\\y=Approx2(zz,r);
y=bestappr(zz,r);
print("r ",r);
print("Approx ",y);
yy=vector(n);
for(i=2,n,yy[i]=y[i-1]);
yy[1]=-(y*AA*y~)/(2*yy[2]);
return(yy);
}

{
Final3(A,B,z,r)=
local(n,AA,zz,y,yy);
n=#A;
AA=vecextract(A,Str(2".."n),Str(2".."n));
zz=vecextract(z,Str(2".."n));
\\y=Approx2(zz,r);
y=bestappr(zz,10^r);
print("r ",r);
yy=vector(n);
for(i=2,n,yy[i]=y[i-1]);
yy[1]=-(y*AA*y~)/(2*yy[2]);
return(yy);
}

{
Final(A,B,z)=
local(n,AA,zz,y,yy);
n=#A;
AA=vecextract(A,Str(2".."n),Str(2".."n));
zz=vecextract(z,Str(2".."n));
\\y=Approx2(zz,r);
\\y=bestappr(zz,10^r);
\\print("r ",r);
\\print("Approx ",y);
yy=vector(n);
for(i=2,n,yy[i]=zz[i-1]);
yy[1]=-(zz*AA*zz~)/(2*yy[2]);
return(yy);
}




