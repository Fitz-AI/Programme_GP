\\ Ce programme calcule une solution commune de A et B, A != 0
Qfn(A,B,r)=
{ my( n = #A, a, b, P, P1, A2, B2, m, mm, cpt, P0, A0, B0, X);

  if( qfsign(A)[1] < floor(n/2),
    [a,b] = dichoqfsign(A,B); \\ on cherche une forme équilibrée
    A = a*A + b*B
  );
  P = Hyperbolique(A);
  A1 = P*A*P~;
  B1 = P*B*P~;
  for(i = 1, n,
    if( A1[i,i] != 0, m=i-1; break)
  );
  \\ b_1,...,b_m sont isotropes, pas b_{m+1}
  P1 = Permute(A1);
  P0 = P = P1*P;
  A0 = A2 = P1*A1*P1~;
  B0 = B2 = P1*B1*P1~;
  mm = m/2;
  \\ On décompose A en plan hyperboliques
  while (1,
    my(C);
    if (cpt > 0, \\ Si on n'a pas de solutions on calcule un autre SETIM
      my(S);
      P = P0;
      A2 = A0;
      B2 = B0;
      S = isotrope2(A2,r);
      P = S*P;
      A2 = S*A2*S~;
      B2 = S*B2*S~;
    );
    C = B2[1..mm, 1..mm];
    if (mm > 5,  mm = 5; C = C[1..5,1..5]);
    if (#C == 5,
      C /= content(C);
      if (type(X = myqfsolve(C,,1)) != "t_INT", break);
    ,
      my(R);
      R = factor(matdet(C), 10^6)[,1];
      if ((#R == 0 || ispseudoprime(R[#R]))
          && type(X = qfsolve(C)) != "t_INT", break);
      print1(".");
    );
    cpt++;
  );    
  print("compteur ", cpt );
  if (type(X) == "t_MAT", X = X[,1]);
  return(X~ * P[1..mm,]);
}
