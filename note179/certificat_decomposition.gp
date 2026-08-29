\\ Certificat de decomposition : verifie EXACTEMENT (arithmetique rationnelle,
\\ pas de flottant) que P = Hyperbolique(Q) realise bien
\\   P * Q * P~ = H^3 (+) R   avec  det(P) != 0,
\\ pour chacun des 12 representants de P^1(Z/8), sur le pinceau du tirage 179.


read("../libre/qfsolve.gp");
read("../programme/Changements_de_bases.gp");

A179=[-6, 0, 2, 10, 1, -3, -6, 8, 8; 0, -3, -2, -9, 8, 8, -5, 1, -2; 2, -2, -2, 8, -2, 0, 1, 2, -8; 10, -9, 8, -6, -6, 2, -2, -4, 3; 1, 8, -2, -6, 10, 5, 4, -1, -7; -3, 8, 0, 2, 5, -4, 2, -4, -9; -6, -5, 1, -2, 4, 2, 4, -3, 7; 8, 1, 2, -4, -1, -4, -3, 10, 3; 8, -2, -8, 3, -7, -9, 7, 3, -8];
B179=[-1, -2, -6, -3, -5, -6, -2, -2, 3; -2, -1, 1, -9, 0, -1, -10, 1, -6; -6, 1, 1, -5, -5, -6, 10, 1, -6; -3, -9, -5, -4, 10, 0, -1, -8, 4; -5, 0, -5, 10, 7, -8, -9, -8, 3; -6, -1, -6, 0, -8, -5, -3, 3, -3; -2, -10, 10, -1, -9, -3, 3, 3, 3; -2, 1, 1, -8, -8, 3, 3, 7, 9; 3, -6, -6, 4, 3, -3, 3, 9, 5];

pts8=[[0,1],[1,1],[2,1],[3,1],[4,1],[5,1],[6,1],[7,1],[1,0],[1,4],[1,2],[3,2]];

\\ verifie que Qd est EXACTEMENT bloc-diagonal H^3 (+) R : diagonale nulle sur
\\ les 6 premieres coordonnees, chaque paire (2k-1,2k) non degeneree

verifie_bloc(Qd)={
  my(n=matsize(Qd)[1], ok=1, detail="");
  for(k=1,3,
    if(Qd[2*k-1,2*k-1]!=0 || Qd[2*k,2*k]!=0,
      ok=0; detail=Str(detail," diag(paire ",k,") non nulle;"));
    if(Qd[2*k-1,2*k]==0,
      ok=0; detail=Str(detail," paire ",k," degeneree (coef nul);"));
  );
  for(i=1,6,
    for(j=1,9,
      my(meme_paire = (ceil(i/2)==ceil(j/2)));
      if(!meme_paire && Qd[i,j]!=0,
        ok=0; detail=Str(detail," entree croisee (",i,",",j,")=",Qd[i,j]," non nulle;"));
    );
  );
  my(R=matrix(3,3,a,b,Qd[6+a,6+b]), R_entier=1);
  for(a=1,3, for(b=1,3, if(type(R[a,b])!="t_INT", R_entier=0)));
  if(!R_entier,
    ok=0; detail=Str(detail," R non entier (denominateur present);"));
  if(R_entier && matdet(R)%2==0,
    ok=0; detail=Str(detail," det(R) pair;"));
  [ok,detail];
}

certificat_pinceau(A,B,nom,pts)={
  printf("--- %s ---\n", nom);
  my(nok=0, nfail=0);
  for(i=1,#pts,
    my(lambda=pts[i][1], mu=pts[i][2], Q=lambda*A+mu*B);
    my(dQ=matdet(Q));
    if(dQ!=0,
      my(P=Hyperbolique(Q));
      my(dP=matdet(P));
      if(dP==0,
        nfail++;
        printf("  (%d:%d) !!! ECHEC : det(P) = 0, P non inversible !!!\n", lambda, mu);
      ,
        my(Qd=P*Q*P~);
        my(chk=verifie_bloc(Qd));
        if(chk[1],
          nok++;
          my(R=matrix(3,3,a,b,Qd[6+a,6+b]));
          printf("  (%d:%d) OK : P*Q*P~ = H^3 (+) R exactement, det(P)=%s, R=%s\n",
                 lambda, mu, dP, R);
        ,
          nfail++;
          printf("  (%d:%d) !!! ECHEC decomposition :%s\n", lambda, mu, chk[2]);
        );
      );
    );
  );
  printf("  %s : %d certificats valides, %d echecs (sur %d points)\n\n", nom, nok, nfail, #pts);
}

print("=== Certificat exact de decomposition P*Q*P~ = H^3 (+) R ===\n");
certificat_pinceau(A179,B179,"tirage 179",pts8);
quit
