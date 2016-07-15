/* Ces algorithmes permettent de calculer un vecteur rationnel x tel que
 * xAx~=0 et xBx~<0 */

DoubleWitt(A,B)=
{ my(n = #A,z,P,B2,Q);

  z = finalsol(A,B);
  P = matid(n);
  P[1,] = z;
  P[2,] /= (P*A*P~)[1,2];
  P = Mordell2(P*A*P~) * P;
  B2 = P*B*P~;
  Q = matid(n);
  Q[3,] /= B3[1,3];
  Q = Mordell6(Q*BB*Q~) * Q;
  return(Q * P);
}

PosNeg(A,B,c)=
{ my(n = #A, A2, B2, d, y, z);

  A2 = A[2..n, 2..n];
  B2 = B[2..n, 2..n];
  d = sign( real(A2[2,2]) );
  y = vector(n-1,i,1);
  y[2] = c*d;

  while (1,
    my (tmp, e);
    tmp = (y*B2*y~) - y[2]*(y*A2*y~);
    e = sign(real(tmp));
    if (e <= 0, break);
    y[2] = 2*y[2];
  );
  z = concat(-(y*AA*y~)/(2*y[1]), y);
  return(z);
}

\\ FIXME: z[2] = 0 + test if z is rational ?
Final(A,B,z,r) =
{ my(n = #A, A2, z2, y);

  A2 = A[2..n, 2..n];
  z2 = z[2..n];
  y = bestappr(z2, 2^r);
  y = concat(-(y*A2*y~)/(2*y[1]), y);
  return (y);
}
