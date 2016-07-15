\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\                               \\
\\  VERSION DU 07 / 02 / 2011    \\
\\                               \\ 
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\contient les algorithmes de minimisation

DEBUG_minp=0;

{minimisation_2(Q,p)=
	local(dd,Ker,dimker,baseker,i,n,Qbk,A3tmp,A3GS,S,CGB1,swap,GS,GG,QfCGB,Qf,M,flagd2,GA,GB);
	\\on suppose que la dimension du noyau est 2
	if(DEBUG_minp>0,print("Minimisation 2 en ",p););
	n=length(Q);
	Ker=kermodp(Q,p);
	\\essai Ker=matsnf(Q,1)[2];
	\\le noyau est engendré par les 2 premiers vecteurs de la base.
	dimker=Ker[1];
	baseker=Ker[2];
	\\essai baseker=Ker;
	A3tmp=matrix(n-2,n-2);
	Qbk=baseker~*Q*baseker;
	\\print("Qbk=",Qbk);
	for(i=1,n-2,
		for(j=1,n-2,
			A3tmp[i,j]=Qbk[i+2,j+2];
		);
	);
	\\print("A3tmp=",A3tmp);
   \\ print("A3tmp[1,1]%p ::",A3tmp[1,1]%p);
    if(A3tmp[1,1]%p==0,
        \\print("\n \nCherie, avant ça plantait ici !!!\n\n");
        GG=matid(n);
        Qf=GG~*baseker~*Q*baseker*GG;
        Qf[,4]*=p;
	    Qf[,5]*=p;
	    Qf[4,]*=p;
	    Qf[5,]*=p;
	    Qf/=p;
	    M=matid(n);
	    M[4,4]=M[5,5]=p;
	    GG=baseker*GG*M;
	    \\print("Return at flag1");
	    return([Qf,GG]);
    );
	GS=gram_schmidt(A3tmp,p);
	\\print("GS=",GS);
	if(type(GS)=="t_COL",
	    \\print("VEC !!!!");
	    \\GS=GS~;
	    CGB1=completebasis(GS);
    	swap=matid(3);
	    swap[1,1]=swap[3,3]=0;
	    swap[1,3]=swap[3,1]=1;
	    CGB1*=swap;
	    CGB=CGB1;
	    GG=matdiagonalblock([matid(2),CGB]);
	    Qf=	GG~*baseker~*Q*baseker*GG;	
	    \\print("Qfavpmodp=",Qf%p);
	    Qf[,4]*=p;
	    Qf[,5]*=p;
	    Qf[4,]*=p;
	    Qf[5,]*=p;
	    Qf/=p;
	    M=matid(n);
	    M[4,4]=M[5,5]=p;
	    GG=baseker*GG*M;
	    \\print("out of min 2 en ",p);
	   \\ print("return at flag 2");
	    return([Qf,GG]);
	);\\fin if
	\\GS=lift(GS);
\\	print("GS=",GS);
   \\ print("A3tmp[1,1]%p ::",A3tmp[1,1]%p);
    /*for(i=1,3,
        dd=delta_diag(GS[2],i);
        if(dd%p==0, print("delta_",i,"=",dd);print("delta_",i,"=",dd%p," %",p););
    );*/
 /*   if(A3tmp[1,1]%p==0,
        print("\n \nCherie, avant ça plantait ici !!!\n\n");
        GG=matid(n);
        Qf=GG~*baseker~*Q*baseker*GG;
        Qf[,4]*=p;
	    Qf[,5]*=p;
	    Qf[4,]*=p;
	    Qf[5,]*=p;
	    Qf/=p;
	    M=matid(n);
	    M[4,4]=M[5,5]=p;
	    GG=baseker*GG*M;
	    return([Qf,GG]);
    );*/
    dd=delta_diag(GS[2],2);
   /* flagd2=0;
    if(dd%p==0,
        GB=matid(n);
        GB[3,4]=GB[4,3]=1;
        GB[3,3]=GB[4,4]=0;
        printp("GB=",GB);
        flagd2=1;
        GA=matid(3);
        GA[1,1]=GA[2,2]=0;
        GA[2,1]=GA[1,2]=1;
        GS[1]=GS[1]*GA;
        GS[2]=GS[2]*GA;
        
        Qf=GG~*baseker~*Q*baseker*GG;
        print("Qf=",Qf);
        Qf[,4]*=p;
	    Qf[,5]*=p;
	    Qf[4,]*=p;
	    Qf[5,]*=p;
	    Qf/=p;
	    M=matid(n);
	    M[4,4]=M[5,5]=p;
	    GG=baseker*GG*M;
	    return([Qf,GG]);
    );*/
    if(denominator(GS[1])%p == 0 || denominator(GS[2])%p==0, print("\n \nCherie, ça va planter !!!\n/!\\ A FAIRE /!\\\n\n"););
	GS[1]=GS[1]%p;
	GS[2]=GS[2]%p;
	\\printp("det GS[1] :",matdet(GS[1]));
	A3GS=GS[2];
	\\print("A3GS=",A3GS);
	S=poll_sch((A3GS[2,2]/A3GS[1,1])%p,(-A3GS[3,3]/A3GS[1,1])%p,p);
	\\print("fin du poll_ch");
	S=concat(S,1)~;
	\\print("S~*A3GS*S %p::",(S~*A3GS*S)%p);
	CGB1=completebasis(S);
	swap=matid(3);
	swap[1,1]=swap[3,3]=0;
	swap[1,3]=swap[3,1]=1;
	CGB1*=swap;
	CGB=GS[1]*CGB1;
	GG=matdiagonalblock([matid(2),CGB]);
	Qf=	GG~*baseker~*Q*baseker*GG;	
	\\print("Qfavpmodp=",Qf%p);
	Qf[,4]*=p;
	Qf[,5]*=p;
	Qf[4,]*=p;
	Qf[5,]*=p;
	Qf/=p;
	M=matid(n);
	M[4,4]=M[5,5]=p;
	GG=baseker*GG*M;
	\\print("out of min 2 en ",p);
	\\print("return at the end");
	return([Qf,GG]);
}

{minimisation_3(Q,p)=
	local(n,Ker,dimker,baseker,GG);
	\\on suppose que la dimension du noyau modulo p est 3
	if(DEBUG_minp>0,print("Minimisation 3 en ",p););
	n=length(Q);
	Ker=kermodp(Q,p);
	dimker=Ker[1];
	baseker=Ker[2];
	GG=baseker;
	Qf=baseker~*Q*baseker;
	M=matid(n);
	M[4,4]=M[5,5]=p;
	GG=GG*M;
	\\print("G=",GG);
	Qf=M~*Qf*M;
	Qf/=p;
	\\print("fin de minim 3");
	return([Qf,GG]);
}\\fin minimisation 3

{is_minimisable(Q)=
	local(SNF,p,fact,nfact,i,Min,Qf,G,n,U,V,SNFF);
    \\print("Qm=",Q);
	SNFF=matsnf(Q,1);
	SNF=SNFF[3];
	U=SNFF[1];
	V=SNFF[2];
	n=length(Q);
	\\printp("SNF=",SNF);
	if(SNF[2,2]==1, return(0););
	if(SNF[3,3]==1,
		\\là, nous sommes dans le cas ou le noyau est de dimension 2 modulo SNF[2].
		\\donc il faut minimiser sur les premiers qui divisent SNF[2] et qui ne sont pas 
		\\dans la base de factorisation.
		snf2=SNF[2,2];
		/*if(is_smooth(snf2),
			print("snf2 est Bf-smooth, pas de minimisation.");
			return(0);
		);*/
		\\print("Minimisation de dimension 2");
		fact=factor(abs(snf2),0);\\pour le moment, on le fait morceau par morceau, factor(*,0) ne teste que les premiers connus par gp
		nfact=length(fact[,1]);
		Qf=Q;
		G=matid(n);
		i=nfact;
		while(i!=0,
		    p=fact[i,1];	
		    if(p!=2 && gcd(p,Pf)==p,i--;next;);
		    Min=minimisation_2(Qf,p);
		    Qf=Min[1];
		    G=G*Min[2];
		    i--;
	    );
	    \\print("out of min 2");
		return([Qf,G]);
	);\\fin if SNF[3,3].
	if(SNF[4,4]==1,
		\\là, nous sommes dans le cas ou le noyau est de dimension 3 modulo SNF[3].
		\\donc il faut minimiser sur les premiers qui divisent SNF[3] et qui ne sont pas 
		\\dans la base de factorisation.
		snf3=SNF[3,3];
		/*if(is_smooth(snf3),
			print("snf3 est Bf-smooth, pas de minimisation.");
			return(0);
		);*/
		/*print("Minimisation de dimension 3");
		snf4=SNF[4,4];
		fact=factor(abs(snf3),0);
		nfact=length(fact[,1]);
		Qf=Q;
		G=matid(n);
		i=nfact;
		while(i!=0,
		    p=fact[i,1];
		    if(gcd(p,Pf)==p,i--;next;);
		    Min=minimisation_3(Qf,p);
		    Qf=Min[1];
		    G=G*Min[2];
		    i--;
	    );
	    return([Qf,G]);*/
		if(DEBUG_minp>0,print("Minimisation 3 en ",snf3););
        G=matid(5);
        G[4,4]=G[5,5]=SNF[3,3];
        G=V*G;
        \\print("G=",G);
        Qf=Q;
        Qf=G~*Q*G;
        Qf/=SNF[3,3];
        \\print("3-Qf=",denominator(Qf));
        return([Qf,G]);
    );\\fin if SNF[3] && SNF[4].
    if(SNF[5,5]==1,
        \\noyau de dimension 4
    	/*if(is_smooth(SNF[4,4]),
			print("snf4 est Bf-smooth, pas de minimisation.");
			return(0);
		);*/
		if(DEBUG_minp>0,print("Minimisation 4 en ",SNF[4,4]););
        G=matid(5);
        G[5,5]=SNF[4,4];
        G=V*G;
        \\print("G=",G);
        Qf=Q;
        Qf=G~*Q*G;
        Qf/=SNF[4,4];
        return([Qf,G]);
        
    ,\\sinon
        \\noyau de dimension 5
        \\on divise toute la matrice par SNF[5,5].
        if(DEBUG_minp>0,print("Minimisation 5 en ",SNF[5,5]););
        G=matid(5);
        G*=SNF[5,5];
        Qf=Q;
        Qf/=SNF[5,5];
        return([Qf,G]);
    );\\fin if SNF[5,5]
    
    
		return([Qf,G]);
	
}

\\minimisation en 2 lorsque le noyau est de dimension 1 ie : reduction de la partie paire
\\on suppose que valuation(det Q,2)>=2
{minim_en_2(Q)=
    local(V,G,n,Qf,detQ,v2,H,Sw,B,Qd,SNFF,SNF);
    if(DEBUG_minp>0,print("Reduction de la partie paire"););
    detQ=matdet(Q);
    v2=valuation(detQ,2);
    n=length(Q);
    G=matid(n);
    Qf=Q;
    G[1,1]*=1/2^(v2\2);
    V=matsnf(Qf,1)[2];
    \\print("V=",V);
    G=V*G;
    Qf=G~*Q*G;
    \\print("Qf=",Qf);
    \\print("G=",G);
    Qd=Qf;
    v2=valuation(matdet(Qf),2);
    if(v2!=0,
        \\print("il en reste :: ", v2);
    ,\\sinon
   if(DEBUG_minp>0, printp("\nstep1A valuation en 2 :",valuation(matdet(Qf),2)););
    return([Qf,G]);        
    );
    \\print(Qf);
    V=matsnf(Qf,1)[2];
    G=G*V;
    Qf=V~*Qf*V;
    \\printp("\nstep1 valuation en 2 : ",valuation(matdet(Qf),2));
    \\on met un 2 en [2,2]
    H=0;
    if(Qf[2,2]%2==1,
        if(Qf[3,3]%2==1,
            H=[0,1,1,0,0]~;
        ,\\sinon, Qf[3,3]%2==0
            H=[0,0,1,0,0]~;
        );
    );
    if(type(H)=="t_COL",
        B=completebasis(H);
        Sw=matrice_swap(2,5);
        B=B*Sw;
        G=G*B;
        Qf=B~*Qf*B;
    );
\\    printp("\nstep2 ",Qf%2,"\nvaluation en 2 :",valuation(matdet(Qf),2));
    Qf=2*Qf;
\\    printp("\nstepmult ",Qf%2,"\nvaluation en 2 :",valuation(matdet(Qf),2));
    V=matid(5);
    V[1,1]=1/2;
    G=G*V;
    Qf=V~*Qf*V;
\\    printp("\nstep3 ",Qf%2,"\nvaluation en 2 :",valuation(matdet(Qf),2));
    V=matid(5);
    V[2,2]=1/2;
    G=G*V;
    Qf=V~*Qf*V;
\\    printp("\nstep4 ",Qf%2,"\nvaluation en 2 :",valuation(matdet(Qf),2));
\\    print("\nQf ",Qf);
    SNFF=matsnf(Qf,1);
    SNF=SNFF[3];
    V=SNFF[2];
    G=G*V;
    Qf=V~*Qf*V;
\\    printp("\nstep5 ",Qf%2,"\nvaluation en 2 :",valuation(matdet(Qf),2));
    if(Qf[1,1]%4==0,
        V=matid(5);
        V[1,1]=1/2;
        G=G*V;
        Qf=V~*Qf*V;
        \\printp("\nstep5fin \nvaluation en 2 :",valuation(matdet(Qf),2));
        return([Qf,G]);
    );
    \\On doit mettre un 4 en haut à gauche
    if((Qf[1,1]/2)%2==1,
        if((Qf[2,2]/2)%2==1,
            H=[1,1,0,0,0]~;
        ,\\sinon, (Qf[2,2]/2)%2==0
            H=[0,1,0,0,0]~;
        );
    );
    B=completebasis(H);
    Sw=matrice_swap(1,5);
    B=B*Sw;
    G=G*B;
    Qf=B~*Qf*B;
\\    printp("\nstep6 ",Qf%4,"\nvaluation en 2 :",valuation(matdet(Qf),2));
    V=matid(5);
    V[1,1]=1/2;
    G=G*V;
    Qf=V~*Qf*V;
    if(DEBUG_minp>0,printp("\nstep7 \nvaluation en 2 :",valuation(matdet(Qf),2)););
    return([Qf,G]);\\On a Qf=2*G~*Q*G
}\\fin minim_en_2


\\procedure qui donne une matrice de swaper
{matrice_swap(i,j)=
    local(M);
    M=matid(5);
    M[i,i]=M[j,j]=0;
    M[i,j]=M[j,i]=1;
    return(M);
}








