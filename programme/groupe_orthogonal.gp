\\ Calcule un élément du groupe orthogonal de A (réflexion)
s(a, x, A)=
{ my(y = a*A, t = 2 * (y*x~) / (y*a~));

  return(x - t*a);
}

\\ cas particulier y = a*A, x = i-ème vecteur base canonique
si(a, i, y) =
{ my(t = 2 * y[i] / (y*a~));

  return(x - t*a);
}

SOn(A,r)=
{ my(n = #A, P = matid(n), S = matrix(n,n), a, y);

  until (y*a~ != 0,
    a = vector(n,i,random([-r,r]));
    y = a*A
  );
  for (i = 1, n, S[i,] = si(a, i, y));
  return(S);
}

\\ Calcule un �l�ment du radical unipotent du groupe orthogonal de A.
\\ A symétrique, équilibrée (|p-m| <= 1, [p,m] = qfsign(A))
isotrope(A, r)=
{ my(n = #A, r2, P=matid(n), m, J, A2, U, E, T, C);

  for(i = 1, n,
    if(A[i,i] != 0, r2=i-1;break)
  );
  m = n - 2*r2;
  J = matrix(r2,r2, i,j, i+j == r2+1);
  A2 = A[r2+1..n-r2, r2+1..n-r2];
  U = matrix(m, r2, i, j, random([-r,r]));
  E = matrix(r2,r2, i, j, random([-r,r]));
  E = E - E~;
  T = -J * U~ * A2;
  C = (T*U)/2 + J*E;
  P = matconcat([matid(r2),0,0; U,matid(m),0; C,T,matid(r2)]);
  return(P~);
}

\\ Calcule un élément d'un sous-groupe parabolique du groupe orthogonal de A

isotrope2(A,r)=
{ my( n = #A, r2, P=matid(n), m, J, J2, A2, F, D, U, E, T, C);

  for(i = 1, n,
    if(A[i,i] != 0, r2=i-1;break)
  );
  m = n - 2*r2;
  J = matrix(r2,r2, i,j, i+j == r2+1);
  A2 = A[r2+1..n-r2, r2+1..n-r2];
  while(1,
    F = matrix(r2,r2,i,j,random([-r,r]));
    if (matdet(F), break)
  );
  J2 = J~*((F~)^-1);
  D = J2*J;
  U = matrix(m, r2, i, j, random([-r,r]));
  E = matrix(r2,r2, i, j, random([-r,r]));
  E = E - E~;
  T = - J2*U~*A2;
  C = (T*U)/2 + J2*E;
  P = matconcat([F,0,0; U,matid(m),0; C,T,D]);
  return(P~);
}
