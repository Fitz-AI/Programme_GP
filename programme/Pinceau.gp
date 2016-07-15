\\ Ce programme renvoie deux entiers u et v tels que uA+vB soit Ã©quilibrÃ©.
dichoqfsign(A,B)=
{ my(n = #A, pol, sol, s, m, s1, t, a, b, c, d, q);

  pol = matdet(A*'x + B);
  s = polrootsreal(pol);
  m = #s;
  s1 = max( abs(floor(s[1])) , abs(ceil(s[m]))) + 1;
  [a,b] = [-s1, s1 + 1];
  c = qfsign(a*A+B)[1];
  d = n - c;
  \\dichotomie
  while (1,
    t = (a+b) / 2;
    d = qfsign(t*A+B)[1];
    if( (ceil(n/2) - c) * (ceil(n/2) - d) < 0, b = t, a = t);
    if (t != 0 && abs(d - n/2) <= 1, break);
  );
  q = denominator(t);
  for (i = 1, oo,
    my(t2 = bestappr(t, 2^i));
    if (denominator(t2) > q, break);
    if (abs(qfsign(t2*A + B)[1]) - n/2 <= 1, t = t2; break);
  );
  return ([numerator(t), denominator(t)]);
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
