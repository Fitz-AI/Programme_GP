\\ Ce programme calcule une solution commune de A et B

Qfn(A,B,r)=
{my( n = #A, a, b, P, P1, A2, B2, m, mm,  d, cpt, Q, A0, B0, C, ZZ, X, Y);

 if( qfsign(A)[1] < floor(n/2),
		  [a,b] = dichoqfsign(A,B); \\ on cherche une forme équilibrée
                   A = a*A+b*B);
  P = Hyperbolique(A);
  A1 = P*A*P~;
  B1 = P*B*P~;
  for(i = 1, n,
            if( A1[i,i] != 0, m=i-1; break)
     );
  P1 = Permute(A1);
  P = P1*P;
  A2 = P1*A1*P1~;
  B2 = P1*B1*P1~;
  mm = m/2;
  Q = P;
  A0 = A2;
  B0 = B2;              \\ On décompose A en plan hyperboliques
  while( d == 0,
                if (cpt > 0,  \\ Si on n'a pas de solutions on calcule un autre SETIM
		            my(S);
		            P=Q;
		            A2 = A0;
		   	    B2 = B0;
		   	    S = isotrope2(A2,r);
		  	    P = S*P;
		  	    A2 = S*A2*S~;
		            B2 = S*B2*S~;
			    );
	 C = B2[1..m,1..m];
         if (mm > 5,  mm = 5; C = C[1..5,1..5]);
         if( #C == 5,
                     C=C/content(C);
                     if( type(X = myqfsolve(C)) != "t_INT",
		                                           d = 1,
		                                          cpt=cpt+1;
	          	) 
           );
         if (mm < 5 ,
	              my(R);
                      R=factor(matdet(C),10^6)[,1];
                      if(matdet(C) == 1  && type(X=qfsolve(C)) != "t_INT", d=1 , cpt=cpt+1;);
                      if(ispseudoprime(R[#R]) == 1  && type(X=qfsolve(C)) != "t_INT",
		                                                                     d = 1,
										     print1(".");
										     cpt = cpt+1;)
             );
	 );    
  print("compteur ", cpt );
  if(type(X) == "t_MAT", X=X[,1]);
  X = X~;
  Y = vector(n);
  for(i = 1, mm, Y[i] = X[i]);
  return(Y*P);
}
