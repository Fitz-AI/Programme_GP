\\ Audit complet : pour les 12 points mod 8 (p=2), sur les trois pinceaux
\\ (temoin, 122, 179), calcule le residu reel de dimension 3 (via Hyperbolique,
\\ ), puis compare DEUX methodes independantes :
\\  (a) mon calcul eps_2 (Prop 10.1.1 + Lemme 10.2.3, avec le correctif p=2) ;
\\  (b) l'arbitre par recherche exhaustive mod 64.


read("../libre/qfsolve.gp");
read("../programme/Changements_de_bases.gp");

countr(M)={my(n=matsize(M)[1], r=0);
  while(2*r+2<=n && M[2*r+1,2*r+1]==0 && M[2*r+2,2*r+2]==0, r++);
  r;}

isotrope_mod64(M)={
  my(trouve=0);
  for(x=0,63,
    if(!trouve,
      for(y=0,63,
        if(!trouve,
          for(z=0,63,
            if(!trouve && (x%2==1 || y%2==1 || z%2==1),
              my(v=[x,y,z]~);
              my(val=(v~*M*v)%64);
              if(val==0, trouve=[x,y,z]);
            );
          );
        );
      );
    );
  );
  trouve;
}

eps2_residu(Q)={
  my(dQ=matdet(Q));
  if(dQ==0, return(-9));
  my(epsQ=iferr(QfWittinvariant(Q,2), E, 0), ok=(epsQ!=0));
  for(tries=1,30,
    if(!ok,
      my(sigma=numtoperm(9,random(362880)));
      my(perm=vecextract(Q,sigma,sigma));
      epsQ=iferr(QfWittinvariant(perm,2), E, 0);
      if(epsQ!=0, ok=1);
    );
  );
  if(!ok, return(-9));
  my(c3=hilbert(-1,dQ,2), eps_q2=epsQ*c3, det_q2=-dQ);
  my(target=hilbert(-1,-det_q2,2));
  if(eps_q2==target, 1, 0);
}

pts8=[[0,1],[1,1],[2,1],[3,1],[4,1],[5,1],[6,1],[7,1],[1,0],[1,4],[1,2],[3,2]];

audit_pinceau(A,B,nom,pts)={
  printf("--- %s ---\n", nom);
  my(naccord=0, ndesaccord=0);
  for(i=1,#pts,
    my(lambda=pts[i][1], mu=pts[i][2], Q=lambda*A+mu*B);
    if(matdet(Q)!=0,
      my(t0=getabstime(), P=0);
      my(Err=alarm(5, P=Hyperbolique(Q)));
      if(type(Err)!="t_ERROR",
        my(Qd=P*Q*P~, r=countr(Qd));
        if(r==3,
          my(resid=matrix(3,3,i2,j2,Qd[2*r+i2,2*r+j2]));
          my(a=eps2_residu(Q));
          my(b=isotrope_mod64(resid));
          my(b_verdict=if(b==0,0,1));
          my(accord=(a==b_verdict));
          if(accord, naccord++, ndesaccord++);
          printf("  (%d:%d) r=3, eps2=%s, mod64=%s -> %s\n", lambda, mu,
                 if(a==1,"isotrope","anisotrope"),
                 if(b==0,"anisotrope",Str("isotrope ",b)),
                 if(accord,"OK","### DESACCORD ###"));
        ,
          printf("  !!! SURPRISE (%d:%d) r=%d != 3 (attendu 3 par la borne generale) !!! a examiner, ne pas ignorer\n", lambda, mu, r);
        );
      ,
        printf("  (%d:%d) timeout Hyperbolique\n", lambda, mu);
      );
    );
  );
  printf("  %s : %d accords, %d desaccords\n\n", nom, naccord, ndesaccord);
}

A179=[-6, 0, 2, 10, 1, -3, -6, 8, 8; 0, -3, -2, -9, 8, 8, -5, 1, -2; 2, -2, -2, 8, -2, 0, 1, 2, -8; 10, -9, 8, -6, -6, 2, -2, -4, 3; 1, 8, -2, -6, 10, 5, 4, -1, -7; -3, 8, 0, 2, 5, -4, 2, -4, -9; -6, -5, 1, -2, 4, 2, 4, -3, 7; 8, 1, 2, -4, -1, -4, -3, 10, 3; 8, -2, -8, 3, -7, -9, 7, 3, -8];
B179=[-1, -2, -6, -3, -5, -6, -2, -2, 3; -2, -1, 1, -9, 0, -1, -10, 1, -6; -6, 1, 1, -5, -5, -6, 10, 1, -6; -3, -9, -5, -4, 10, 0, -1, -8, 4; -5, 0, -5, 10, 7, -8, -9, -8, 3; -6, -1, -6, 0, -8, -5, -3, 3, -3; -2, -10, 10, -1, -9, -3, 3, 3, 3; -2, 1, 1, -8, -8, 3, 3, 7, 9; 3, -6, -6, 4, 3, -3, 3, 9, 5];

A122=[9, -4, -10, -2, -2, 8, -5, 3, 3; -4, 6, -7, 0, 1, -10, 6, -4, 0; -10, -7, -6, -5, 7, 8, 8, 6, -7; -2, 0, -5, -6, 0, -10, -5, 1, -2; -2, 1, 7, 0, 1, -8, -1, -4, 8; 8, -10, 8, -10, -8, -1, -1, -1, 10; -5, 6, 8, -5, -1, -1, -1, 3, -2; 3, -4, 6, 1, -4, -1, 3, 8, -5; 3, 0, -7, -2, 8, 10, -2, -5, -10];
B122=[-2, -5, 7, -6, 10, 7, -3, 5, 8; -5, -9, 10, -5, -5, 6, -6, -7, -2; 7, 10, -2, 8, 4, -9, -1, 7, 2; -6, -5, 8, -4, -2, 8, 10, 6, 6; 10, -5, 4, -2, 8, 1, -6, 0, 0; 7, 6, -9, 8, 1, 1, 2, 0, 6; -3, -6, -1, 10, -6, 2, 4, -5, -7; 5, -7, 7, 6, 0, 0, -5, 6, -6; 8, -2, 2, 6, 0, 6, -7, -6, 3];

symrand(n,H)={my(M=matrix(n,n));for(i=1,n,for(j=i,n,my(c=random(2*H+1)-H);M[i,j]=c;M[j,i]=c));M;}
setrand(999001);
Avalid=symrand(9,10); Bvalid=symrand(9,10);

print("=== Audit croise eps_2 vs arbitre mod 64, granularite mod 8 ===\n");
audit_pinceau(Avalid,Bvalid,"pinceau temoin",pts8);
audit_pinceau(A122,B122,"tirage 122 (temoin connu (41,61))",pts8);
audit_pinceau(A179,B179,"tirage 179 (suspect)",pts8);
quit
