\\ Arbitre: isotropie d'une forme ternaire sur Q_2 par
\\ recherche exhaustive d'un vecteur primitif v (pas tous pairs) mod 2^6=64
\\ tel que v~*M*v = 0 mod 64. 

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

\\ ===== Validation sur trois cas connus =====

print("=== Validation de l'arbitre mod 64 ===");
M1=[1,0,0;0,1,0;0,0,1];
print("x^2+y^2+z^2 (doit etre ANISOTROPE) : ", isotrope_mod64(M1));
M2=[1,0,0;0,1,0;0,0,-2];
print("x^2+y^2-2z^2 (doit etre ISOTROPE, ex. (1,1,1)) : ", isotrope_mod64(M2));
\\ x^2+2xy+3y^2-z^2 : non diagonal, det=-1*(3-1)
\\ isotrope sur Q (donc sur Q_2) via x=1,y=0,z=1
M3=[1,1,0;1,3,0;0,0,-1];
print("x^2+2xy+3y^2-z^2, non diagonal (doit etre ISOTROPE, ex. (1,0,1)) : ", isotrope_mod64(M3));
quit
