
Alea(n,m)=
{ my(A);

  A = matrix(n, n, i, j, random([-m,m]));
  A=A+A~;
  A=A\2;
  return(A);
}
