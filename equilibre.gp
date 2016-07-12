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
