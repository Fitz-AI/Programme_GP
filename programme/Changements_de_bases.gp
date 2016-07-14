/*  Les prochains programmes permettent simplement de faire des changements de variables pour décomposer
 une forme quadratique en plans hyperboliques. Il est possible de trouver plus de précision dans 
mon article Effective Hasse principle for the intersection of two quadrics. */

completebasis(v,redflag=0) =
{ my(U,n,re);

  U = (mathnf(Mat(v~),1)[2]~)^-1;
  n = length(v~);
  if( n==1 || !redflag, return(U));
  re = qflll(vecextract(U,1<<n-1,1<<(n-1)-1));
  return( U*matdiagonalblock([re,Mat(1)]));
}

Mordell2(A)=
{ my(n = #A, P = matid(n));

  for(i = 3, n, P[i,2] = -A[i,1]/A[2,1]);
  return(P);
}


Mordell3(A)=
{ my(n = #A, P = matid(n), S);

  for(i = 2, n, P[i,1]=-A[i,2]);
  A2 = P*A*P~;
  S = matid(n);
  S[2,] = 2*S[2,]-A2[2,2]*S[1,];
  S[1,1] = 1/2;
  P = S*P;
return(P);
}

Mordell6(A)=
{ my(n = #A, P = matid(n));

  P[2,3] = -A[2,1];
  for(i = 4, n, P[i,3] = -A[i,1]);
  return(P);
}  


BaseWitt(A)=
{ 
my(n = #A, P=matid(n), X, tmp);

  X=qfsolve(A);
  if(type(X)=="t_MAT", X=X[,1]);
  if(type(X)=="t_INT",return(1));
  P=completebasis(X)~;
  tmp=P[1,];P[1,]=P[n,];P[n,]=tmp;
  return(P);
}


Redqf(A)=
{ my(P, A2, P2, g);

  P = BaseWitt(A);
  if (type(P) == "t_INT", return(1));
  A2 = P*A*P~; 
  g = A2[1,2];
  A2[1,] /= g;
  A2[,1] /= g;
  P[1,] /= g;
  P2 = Mordell2(A2);
  A2 = P2*A2*P2~ ;
  P = P2*P;
  P2 = Mordell3(A2);
  P=P2*P;
  return(P);
}  

 
Hyperbolique(A)=

{ my( n = #A, r, s, u, i=1, t=-1, P=matid(n));

  [r,s] = qfsign(A);
  u= min(r,s);
  while(t<u,
	    my(A2,tmp,C=matid(n));
	    A2 = A[i..n,i..n];
            tmp = Redqf(A2);
            if(type(tmp) == "t_INT", break);
            for(j = i, n,
		  for(k = i, n, C[j,k] = tmp[j-i+1,k-i+1])
		);
            P = C*P; 
            A = C*A*C~;
            i += 2;
  	    t += 1;
    );
  return(P);
}


Permute(A)=
{ my(n = #A, P=matrix(n,n), j, k=1);

  P[1,1]=1;
  for ( i = 2, n, 
               if ( i % 2 == 0,
		               P[i,n-j] = 1 ;
		               j = j+1 ,
                               P[i,1+k] = 1 ;
		               k = k+1
		  );
      );
  return(P~);
}  
