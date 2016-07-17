\\\\\\\\\\\\\\\\\\\\\\\\\\
\\
\\
\\
\\
\\
\\
\\\\\\\\\\\\\\\\\\\

\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
\\                               \\
\\  VERSION DU 31 / 01 / 2011    \\
\\                               \\ 
\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


\\procédé de gram-schmidt
{gram_schmidt(Q,p)=
	local(n,G,u,v,i,nu,nv,res,M,V,val);
	n=length(Q);
	M=matid(n);
	G=matrix(n,n);
	u=M[,1];
	nu=u~*Q*u;
	\\print(nu);
	G[,1]=u;
	for(i=2,n,
		\\print("G=",G);
		\\print("");
		u=M[,i];
		for(j=1,i-1,
		    val=G[,j]~*Q*G[,j];
		    if(gcd(denominator(val),p)!=1,print("zop"));
		    if(val%p==0,
		        \\print("\nVecteur orth !!!!!\n");
		        \\print("Q=",Q);
		        V=G[,j];
		        V=V%p;
		        return(V)
	        );
	    );
		v=u-sum(j=1,i-1,((G[,j]~*Q*u)/(G[,j]~*Q*G[,j]))*G[,j]);
	    \\v=u-sum(j=1,i-1,(G[,j]~*Q*u)*G[,j]);
		\\v=v/(v~*Q*v);
		\\print("v=",v);
		\\print("");
		G[,i]=v;
	);
	
	\\res=mattranspose(G)*Q*G;
	res=G~*Q*G;
	\\printp(res);
	\\print("G=",G);
	\\print("res=",res);
	return([G,res]);
}

\\donner le delta_i
{delta_diag(Qdiag,i)=
	local(n,j,delta);
	n=length(Qdiag);
	if(i==1,return(Qdiag[1,1]););
	delta=1;
	for(j=1,i,
		\\print("i:", delta);
		delta*=Qdiag[j,j];
	);
	\\print(delta);
	return(delta);
}

addhelp(gram_schmidt,"gram_schmidt(Q,p) : orthogonalisation de Gram-Schmidt de la forme quadratique Q modulo p. Donne la forme reduite Qr et le changement de base B sous la forme [B,Qr]. Si trouve un vecteur de norme nulle modulo p, le renvoie.");
addhelp(delta_diag,"delta_diag(Qdiag,i) : renvoie le produit pour j de 1 a i des coefficients diagonaux de la forme Qdiag.")

