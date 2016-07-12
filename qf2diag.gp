{
qf2diag(A,B)=
local(n,pol,m,sol,boo,K,P);
n=#A;
kill(x);
pol=matdet(A*x+B);
m=polsturm(pol);
sol=polroots(pol);
\\boo=vector(n);
\\for(i=1,n, if (real(sol[i])==sol[i], boo[i]=1; sol[i]=real(sol[i])));
for(i=1,m, sol[i]=real(sol[i]));
K=vector(n);
for(i=1,n, K[i]=matker(A*sol[i]+B)~[1,]);
i=1;
while(i<n,
            \\if ( boo[i]==1 , i=i+1);
            \\if ( boo[i]==0 , K[i+1]=imag(K[i]); K[i]= real(K[i]);  i=i+2);
            if ( i <= m, i=i+1);
            if ( i>m , K[i+1]=imag(K[i]); K[i]= real(K[i]);  i=i+2);
     );
P=matid(n);
for(i=1,n,P[i,]=K[i]);
return([P,m]);
}

{
qf2red(A,B,m)=
local(n,i,P,C,Q);
n=#A;
if (m != n ,i=m+1);
P=matid(n);
if (m != n, while( i <= n,
               C=vecextract(A,Str(i".."i+1),Str(i".."i+1));
               Q=qfjacobi(C)[2]~;
               P[i,i]=Q[1,1]; P[i,i+1]=Q[1,2]; P[i+1,i]=Q[2,1]; P[i+1,i+1]=Q[2,2];
               i=i+2)
   );
return(P);
}


{
finaldiag(A,B)=
local(n,P,m,A2,B2,Q,A3,B3,A4,P2,A5);
n=#A;
[P,m]=qf2diag(A,B);
A2=P*A*P~;
B2=P*B*P~;
Q=qf2red(A2,B2,m);
A3=Q*A2*Q~;
B3=Q*B2*Q~;
P=Q*P;
for(i=1,n, P[i,]=P[i,]/sqrt(abs(A3[i,i])));
A4=P*A*P~;
\\if (m==n,
\\         A5=matid(n);
\\         for(i=1,n,A5[i,i]=A4[i,i]);
\\         P2=qflllgram(A5)~;
\\         P=P2*P;
\\);
return(P);
}

{
solmix(A2,B2)=
local(pol,pol2,pol3,y1,x1,x2,x3,sol);
x='x; y='y; z='z;
\\P=finaldiag(A,B);
\\A2=P*A*P~;
\\B2=P*B*P~;
if(A2[1,1]<0 , A2=-A2);
pol=[x,y,z]*B2*[x,y,z]~;
pol2=substpol(pol,x^2,-sign(A2[2,2])*y^2-sign(A2[3,3])*z^2);
pol3=substvec(pol2,[x,z],[0,1]);
sol=polroots(pol3);
x1= sqrt(-sign(A2[2,2])*sol[1]^2-sign(A2[3,3])*1);
x2= sqrt(-sign(A2[2,2])*sol[2]^2-sign(A2[3,3])*1);
if( real(x1)==x1 , x3=x1 ; y1=sol[1] , x3=x2 ; y1=sol[2]);
\\y1=-B2[2,3]-sign(B2[2,3])*sqrt(B2[2,3]^2-(-sign(A[2,2])*B2[1,1]+B2[2,2])^2);
\\z1=-sign(A[2,2])*B2[1,1]+B2[2,2];
\\x1= sqrt(-sign(A2[2,2])*y1^2-sign(A2[3,3])*z1^2);
return(real([x3,y1,1])); \\ Plus de P
}

{
solcomp(A2,B2)=
local(pol,pol2,pol3,pol4,y1,z1,x1,x,y,z);
\\P=finaldiag(A,B);
\\A2=P*A*P~;
\\B2=P*B*P~;
x='x;y='y;z='z;w='w;
pol=[x,y,z,w]*B2*[x,y,z,w]~;
pol2=subst(pol,x,y);
if ( sign(B2[1,2])==sign(B2[3,4]), pol3=subst(pol2,z,-w) , pol3=subst(pol2,z,w) );
pol4=substvec(pol3,[x,z],[0,0]);
y1=polroots(subst(pol4,w,1))[1];
if ( sign(B2[1,2])==sign(B2[3,4]), z1=-1 , z1=1 );
x1=y1;
return(real([x1,y1,z1,1])); \\ Plus de P
}

{
solreel(A2,B2)=
local(n,tmp1,tmp2,pos,neg,M,m,in,sol);
n=#A2;
x='x;y='y;z='z;
\\P=finaldiag(A,B);
\\A2=P*A*P~;
\\B2=P*B*P~;
for(i=1,n, if ( A2[i,i]>0 , tmp1=tmp1+1));
for(i=1,n, if ( A2[i,i]<0 , tmp2=tmp2+1));
pos=vector(tmp1);
neg=vector(tmp2);
k=0;l=0;
for(i=1,n, if ( A2[i,i]>0 , pos[k+1]=i; k=k+1)); 
for(i=1,n, if ( A2[i,i]<0 , tmp2=tmp2+1; neg[l+1]=i;l=l+1));
M=vecmax(vector(k,i,B2[pos[i],pos[i]]),&v1);
m=vecmin(vector(k,i,B2[pos[i],pos[i]]),&v2);
for(i=1,l, if (-B2[neg[i],neg[i]]< M && -B2[neg[i],neg[i]]> m, v3=neg[i]));
\\if ( v3==0 , for(i=1,l, if (-BB[neg[i],neg[i]]> m, v3=neg[i])));
in=[pos[v1],pos[v2],v3];
A3=matdiagonal([A2[in[1],in[1]],A2[in[2],in[2]],A2[in[3],in[3]]]);
B3=matdiagonal([B2[in[1],in[1]],B2[in[2],in[2]],B2[in[3],in[3]]]);
pol=[x,y,z]*B3*[x,y,z]~;
pol2=substpol(pol,z^2,y^2+x^2);
pol3=subst(pol2,y,1);
[x1,x2]=polroots(pol3);
z1=sqrt(x1^2+1);
z2=sqrt(x2^2+1);
if(real(z1)==z1, z3=z1; x3=x1, z3=z2;x3=x2);
sol=vector(n);
sol[in[1]]=real(x3);sol[in[2]]=1;sol[in[3]]=real(z3);
return(sol); \\ Plus de P
}

{
finalsol(A,B)=
local(n,P,n1,sol,sol2,A2,B2,A3,B3);
n=#A;
n1=polsturm(matdet(A*x+B));
P=finaldiag(A,B); \\ nouveau
A2=P*A*P~;       \\
B2=P*B*P~;        \\
if (n1==n , sol= solreel(A2,B2);return(sol*P));
if (n1==0 , 
           sol=vector(n) ; 
           A3=vecextract(A2,"1..4","1..4");
           B3=vecextract(B2,"1..4","1..4");
           sol2= solcomp(A3,B3);
           for(i=1,4,sol[i]=sol2[i]);
           return(sol*P);
   );
A3=vecextract(A2,Str(n1".."n1+2),Str(n1".."n1+2));
B3=vecextract(B2,Str(n1".."n1+2),Str(n1".."n1+2));
sol=vector(n);
sol2=solmix(A3,B3);
for (i=0,2, sol[n1+i]=sol2[i+1]);
return(sol*P);
}







              



