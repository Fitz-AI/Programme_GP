\\ Ce programme renvoie deux entiers u et v tels que uA+vB soit Ã©quilibrÃ©.
dichoqfsign(A,B)=
{ my(n = #A, pol, sol, s, m, s1, a, b, q, t=1, c, d, u, v);
  \\ test de m=1
  if( qfsign(A+B)[1] >= floor(n/2), return(1));
  \\ test de -1
  if( qfsign(-A+B)[1] >= floor(n/2), return(-1));
  pol = matdet(A*x+B);
  s = polrootsreal(pol);
  m = #s;
  s1 = max( abs(floor(s[1])) , abs(ceil(s[m]))) + 1;
  [a,b] = [-s1,s1];
  \\dichotomie
  while( t==0 || q != ceil(n/2),   
    c = qfsign(a*A+B)[1];
    d = qfsign(t*A+B)[1];
    if( (ceil(n/2) - c) * (ceil(n/2) - d) < 0, b = t, a = t);
    t = (a+b)/2;
    q = qfsign(t*A+B)[1]
  );
  u = numerator(t);
  v= denominator(t);
  return([u,v]);
}

WittPinceau(A,B,r)=
{ my(n = #A, pol, s, m, s1);
  pol = matdet(A*x+B);
  s = polrootsreal(pol);
  m = #s;
  s1 = max( abs(floor(s[1])) , abs(ceil(s[m])));
  while (1, 
    my(a, b, A2, P);
    printf(".");
    a = random([-2^r,2^r]);
    b = random([-2^r,2^r]);
    \\  if ( a == 0, a += 1 );
    \\  b = random( [-abs(a)/s1, abs(a)/s1]); \\ On calcule a et b tels que abs(a/b) > abs(s1)
    A2 = a*A + b*B;
    E = alarm(3, P = Hyperbolique(A2));
    if (type(E) == "t_ERROR", next); 
    A2 = P*A2*P~; \\ On décompose A2 en plans hyperboliques
    if (A2[n-2,n-2] == 0, return([a,b]))
  );
}
