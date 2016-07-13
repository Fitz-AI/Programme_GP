\\ Q[...]
Mordell2(Q)=
{ my(n = #Q, P = matid(n));

  for(i = 3, n, P[i,2] = -Q[i,1]/Q[2,1]);
  return(P);
}

Mordell6(Q)=
{ my(n = #Q, P = matid(n));

  P[2,3] = -Q[2,1];
  for(i = 4,n, P[i,3] = -Q[i,1]);
  return(P);
}


DoubleWitt(A,B)=
{ my(n = #A,z,P,BB,Q);

  z = finalsol(A,B);
  P = matid(n);
  P[1,] = z;
  P[2,] /= (P*A*P~)[1,2];
  P = Mordell2(P*A*P~) * P;

  BB = P*B*P~;
  Q = matid(n);
  Q[3,] = Q[3,]/BB[1,3];
  Q = Mordell6(Q*BB*Q~) * Q;
  return(Q * P);
}

PosNeg(A,B,c)=
{ my(n = #A, AA, BB, d, y, z);

  AA = A[2..n, 2..n];
  BB = B[2..n, 2..n];
  d = sign( real(AA[2,2]) );
  y = vector(n-1,i,1);
  y[2] = c*d;

  while (1,
    my (tmp, e);
    tmp = (y*BB*y~) - y[2]*(y*AA*y~);
    e = sign(real(tmp));
    if (e <= 0, break);
    y[2] = 2*y[2];
  );
  z = concat(-(y*AA*y~)/(2*y[1]), y);
  return(z);
}

\\ FIXME: z[2] = 0 + test if z is rational ?
Final(A,B,z,r) =
{ my(n = #A, AA, zz, y);

  AA = A[2..n, 2..n];
  zz = z[2..n];
  y = bestappr(zz, 10^r);
  print("r = ", r);
  y = concat(-(y*AA*y~)/(2*y[2]), y);
  return (y);
}
