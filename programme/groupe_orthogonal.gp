\\ Calcul un élément du groupe orthogonal de A

s(a,x,A)=
{ my(r);

  r = x-(2 * (a*A*x~) / (a*A*a~)) * a;
  return(r);
}

SOn(A,r)=
{ my(n = #A, P = matid(n), S = matrix(n,n), a = vector(n));

  while ( (a*A*a~)==0 , a=vector(n,i,random([-r,r])));
  for (i = 1, n, S[i,] = s(a, P[i,], A));
  return(S);
}

\\ Calcule un élément du radical unipotent du groupe orthogonal de A 

isotrope(A,r)=
{ my( n = #A, r2, P=matid(n), m, J, A2, U, E, T, C);

  for(i = 1, n,
              if(A[i,i] != 0, r2=i-1;break)
   );
  m = 2*r2;
  J = matrix(r2,r2);
  for (i = 1, r2,
                for(j = 1, r2,
		            if (i+j == r2+1, J[i,j] = 1)
	           )
      );
  A2 = A[r2+1..r2+(n-m),r2+1..r2+(n-m)];
  U = matrix(n-m, r2, i, j, random([-r,r]));
  for (i = r2+1 , n-r2 ,
                     for (j=1 , r2 ,
		                   P[i,j]=U[i-r2,j]
		         )
      );
  E = matrix(r2, r2, i, j, random([-r,r]));
  E = E-E~;
  T = -J*U~*A2;
  for (i = n-r2+1, n ,
                      for (j = r2+1, n-r2,
		                          P[i,j] = T[i-(n-r2),j-r2]
			  )
      );
  C = (-1/2)*J*(U~*A2*U)+J*E;
  for (i = n-r2+1, n,
                   for (j = 1 ,r2, P[i,j]=C[i-(n-r2),j]
		       )
      );
  return(P~);
}

\\ Calcule un élément d'un sous-groupe parabolique du groupe orthogonal de A

isotrope2(A,r)=
{ my( n = #A, r2, P=matid(n), m, J, A2, tmp, F, D, M0, U, E, T, C);

  for(i = 1, n,
              if(A[i,i] != 0, r2=i-1;break)
   );
  m = 2*r2;
  J = matrix(r2,r2);
  for (i = 1, r2,
                for(j = 1, r2,
		            if (i+j == r2+1, J[i,j] = 1)
	           )
      );
  A2 = A[r2+1..r2+(n-m),r2+1..r2+(n-m)];
  while(tmp == 0,
                 F=matrix(r2,r2,i,j,random([-r,r]));
		 tmp=matdet(F);
       );
  for (i = 1 , r2 ,
                  for (j = 1 , r2 ,
		                  P[i,j]=F[i,j];
		      )
      );
  D = J~*((F~)^-1)*J;
  for (i = n-r2+1 ,n ,
                     for (j = n-r2+1, n,
		                       P[i,j]=D[i-n+r2,j-n+r2];
		         )
      );
  M0 = matrix(n-m, n-m, i, j, random([-r,r]));
  for (i = r2+1, n-r2,
                     for (j = r2+1, n-r2,
		                         P[i,j]=M0[i-r2,j-r2];
			 )
      );
  U = matrix(n-m, r2, i, j, random([-r,r]));
  for (i = r2+1, n-r2,
                     for (j = 1, r2, P[i,j]=U[i-r2,j])
      );
  E = matrix(r2, r2, i, j, random([-r,r]));
  E = E-E~;
  T = -J~*((F~)^-1)*U~*A2*M0;
  for (i = n-r2+1, n,
                    for (j = r2+1, n-r2,
		                       P[i,j]=T[i-(n-r2),j-r2]
		        )
      );
  C = (-1/2)*J*((F~)^-1)*(U~*A2*U)+J*((F~)^-1)*E;
  for (i = n-r2+1, n,
                     for (j = 1, r2,
		                    P[i,j] = C[i-(n-r2),j]
		     )
      );
  return(P~);
}