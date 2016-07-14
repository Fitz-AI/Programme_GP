
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\                               \\
\\  VERSION DU 31 / 01 / 2011    \\
\\                               \\ 
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

print("");
print("===== poll_sch.gp =====");
print("- poll_sch(k,m,n)");
print("");

global(prof);
prof=0;
{poll_sch(k,m,n)=
	local(m0,x0,u,v,temp,count,x1,y1,x2,y2,sol2,X,Y,Matr,i,j);
	prof++;
	if(gcd(k,n)!=1 || gcd(m,n)!=1, error("pgcd !"););
\\	print("entree dans poll_ch avec");
\\	print("k=",k," m=",m," n=",n);
	if(abs(n)<20,   
	    for(i=0,abs(n)-1,
	        for(j=0,abs(n)-1,
	            if(((i^2+k*j^2)%n)==m,
                    return([i,j]);
                );
            );
        );
    );
	\\on regarde si il y a une solution triviale
	s=sol_triv(k,m,n);
	if(type(s)=="t_VEC", 
	    \\print("sol triv :",s);
	    return(s);
    );
	sol2=0;
	
	while(sol2==0,
    	\\tirage du m0
	    temp=tirage_m0(k,m,n);
	    m0=temp[1];
	    x0=temp[2];
	    u=temp[3];
	    v=temp[4];
	
	
	
	    \\on récupère une solution de l'équation
	    \\x^2+k*y^2 = m0 mod n
	    \\on la note [x1,y1]
	    sol2 = resoud2(k,m0,n,x0);
    );
	x2=sol2[1];
	y2=sol2[2];
	
	Matr=[x2,y2;u,v];
	\\calcul de la solution de fin
	temp=prod_Qk(Matr,k);
	x1=temp[1];
	y1=temp[2];
	
/*	print("gcd(m0,n)=",gcd(m0,n));
	print("m0=",m0);
	print("gcd(m,n)=",gcd(m,n));
	print("gcd(x1,n)=",gcd(x1,n));
	print("gcd(y1,n)=",gcd(y1,n));*/
	X=lift(Mod(m*x1/m0,n));
	Y=lift(Mod(m*y1/m0,n));
/*	print("sortie de poll_ch avec ",[X,Y]);
	print("les valeurs étaient: k=",k," m=",m," n=",n);
	print("verification : m (mod n):",Mod(m,n));
	print("             x^2+k*y^2  :",Mod(X^2+k*Y^2,n));
	print("");*/
	return([X,Y]);
	
}\\fin poll_sch


{tirage_m0(k,m,n)=
	local(u,v,m0,x0);
	\\print("tirage de m0");
	while(1,
		u=random(abs(n));
		v=random(abs(n));
		\\print("u=",u);
		m0 = (m*(u^2+k*v^2));
		if(m0==0,next;);
		m0=m0%n;
		\\print("m0=",m0);
        if(isprime(m0) && kronecker(-k,m0)==1 && gcd(m0,n)==1,\\ispseudoprime
        	x0=lift(sqrt(Mod(-k,m0)));
        	if(x0^2+k==0,
        		x0+=m0;
            );
            break;
            \\print("Fail !");
        );\\fin if
    );\\fin while
    \\x0=lift(sqrt(Mod(-k,m0)));
    \\print("m0=",m0," x0=",x0," u=",u," v=",v);
    \\print("sortie du tirage de m0");
    return([m0,x0,u,v]);
}


{poll_sch_reduc(k,x0,m0,n)=
	local(xi,mi,i,l,M,h,s,T,M1,U,V);
	\\print("---------------Go for a reduc !");
	l=sqrt(abs(k)); 
	xi=vector(1);
    xi[1]=x0;
    mi=vector(1);
    mi[1]=m0;
	i=1;
	\\print("n=",n," k=",k," m0=",m0," x0=",x0);
    if(k>0,
        while(1,
        	\\print("");
        	\\print("k>0");
        	\\print("mi=",mi);
        	\\print("xi=",xi);
            mi=concat(mi,(xi[i]^2+k)/mi[i]);
\\            print("A mi=",mi);
            if(mi[i+1]==0,/*print("A GAFFE !");*/i--;break;);
            xi=concat(xi,min(lift(Mod(xi[i],mi[i+1])),lift(mi[i+1]-Mod(xi[i],mi[i+1]))));
\\        	print("xi=",xi);
            if(xi[i]<=mi[i+1] && mi[i+1]<=mi[i],
                break;
            );
            i++;
        );\\fin while
    ,\\ si k<0
        while(1,
        	\\print("");
	        \\print("k<0");
        	\\print("mi=",mi);
        	\\print("xi=",xi);
            mi=concat(mi,(xi[i]^2+k)/mi[i]);
\\           \\ print("B mi=",mi);
            if(mi[i+1]==0,/*print("GAFFE !");*/i--;break;);
            xi=concat(xi,min(lift(Mod(xi[i],mi[i+1])),lift(mi[i+1]-Mod(xi[i],mi[i+1]))));
        	\\print("xi%n=",xi%n);
            if(abs(mi[i+1])<=l,
                break;
            );
            i++;
        );\\fin while
    );\\fin if sgn k
\\    print("Fin des calculs de réduction");
    i+=1;\\pour avoir la taille de xi et mi
    M=prod(j=2,i-1,mi[j]^2)*mi[1]*mi[i];
	pol=x^2+k;
	h=prod(j=1,i-1,xi[j]+x);
	h=lift(Mod(h,pol));
	\\print("h=",h);
    s=subst(h,x,0);
    \\print("s=",s);
    t=subst(h,x,1)-s;
    M1=prod(j=2,i,mi[j]);
\\    print("M1=",M1);
\\    print("gcd(M1,n)=",gcd(M1,n));
\\    print("Mod(M1,n)=",Mod(M1,n));
    if(M1==n,print("ouf");return(0););
    if(Mod(M1,n)==0,print("\n/!\\ A FAIRE /!\\"););
  	U=lift(Mod(s/M1,n));
    V=lift(Mod(t/M1,n));
    \\print("U=",U);
    \\print("V=",V);
    \\print("U^2+k*V^2 mod n = ",Mod(U^2+k*V^2,n));
    \\print(" m0/mI mod n = ",Mod(m0/mi[i],n));
    \\print("");
    \\print("mI=",mi[i]," U=",U," V=",V); 
\\    print("------------Fin de la reduc");
	return([mi[i],U,V]);
}\\fin poll_sch_reduc




\\A est un vecteur de taille u,2
\\qui représente les éléments du type
\\A[i]+sqrt(-k)*B[i]
\\Cette procédure fait le produit des éléments et donne le résultat
{prod_Qk(A,k)=
	local(pol,temp,u,v,j);
	\\print("un prod_Qk");
	pol=x^2+k;
	temp=prod(j=1,length(A[,1]),A[j,1]+A[j,2]*x);
	temp=lift(Mod(temp,pol));
	u=subst(temp,x,0);
	v=subst(temp,x,1)-u;
	\\print("fin du prod_Qk");
	return([u,v]);
}

\\une procédure pour voir si il y a des solutions évidentes
{sol_triv(k,m,n)=
	local();
	/*if(m==1,
	    return([1,0]);
    );*/
	if(k==-1,
		if(m==1,
		\\	print("Sol triv:[1,0]");
			return([1,0]);
		);
		if(m==-1,
			\\print("Sol triv:[0,1]");
			return([0,1]);
		);
	);
	if(k==1 && m==1,
	\\	print("Sol triv:[1,0]");
		return([0,1]);
	);
	if(k==3 && m==4 && n==5,
		return([1,1]);
	);
	return(0);
}\\fin sol_triv


\\resoud2
\\on cherche a résoudre l'équation
\\x^2+k*y^2=m0 mod n
{resoud2(k,m0,n,x0)=
	local(temp,mI,U,V,x2,y2,x3,y3,Matr,res);
	\\on commence par effectuer une reduction
	temp=poll_sch_reduc(k,x0,m0,n);
	if(temp==0,return(0));
	mI=temp[1];
	U=temp[2];
	V=temp[3];
	
	if(mI==1,
	\\	print("mI=1");
		return([U,V]);
	);
\\	print("mI!=1");
	\\il nous faut une solution de l'équation:
	\\ x^2 + k*y^2=mI mod n
	\\mais en fait, on va résoudre
	\\x^2 - mI y^2= -k mod n
	y2=0;
	while(y2==0,
    	temp=poll_sch(-mI,-k,n);
	    \\print("Solution avec mI:");
	    \\print("k=",-mI," m=",-k," n=",n);
	    x2=temp[1];
	    y2=temp[2];
	    \\print("x=",x2," y=",y2);
    );
\\	print("gcd(y2,n):",gcd(y2,n));
	\\print("y2=",y2);
	\\on remonte d'un cran
	\\"GERER LA REMONTEE LORSQUE Y==0"
	if(y2==0,print("\n/!\\ A FAIRE /!\\ "););
	/*if(y2==0,
	    temp=poll_sch(-1,-mI,n);
    	x4=temp[1];
    	y4=temp[2];
	    x3=x4;
	    y3=lift(Mod(y4/x2,n));
    ,\\sinon*/
    	y3=lift(Mod(1/y2,n));
    	x3=lift(Mod(x2/y2,n));
	\\);
	
	
	Matr=[x3,y3;U,V];
	res=prod_Qk(Matr,k);
	\\print("Sortie de resoud2");
	return(res);
}


{testpoll(k,m,n)=
	local(sol,res);
	print("");
	sol=poll_sch(k,m,n);
	res=Mod(sol[1]^2+k*sol[2]^2,n);
	print("sol=",sol);
	print("Valeur de m: ",Mod(m,n));
	print("Verification: ",res);
	if(Mod(m,n)!=res,error("Mauvaise solution !"));
	return(sol);
}

addhelp(poll_sch,"poll_sch(k,m,n): Etant donnés k,m et n entiers tels que gcd(k,n)=gcd(m,n)=1, poll_sch(k,m,n) donne [x,y] tel que x^2 + k*y^2 = m (mod n)");


