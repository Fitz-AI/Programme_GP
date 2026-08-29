units=[1,3,5,7];

classe_de(p)={
  my(reps=[]);
  for(i=1,4, my(u=units[i]); reps=concat(reps,[[(u*p[1])%8,(u*p[2])%8]]));
  vecsort(reps);
}

prims=[];
ajoute_si_primitif(a,b)={
  if(a%2==1 || b%2==1, prims=concat(prims,[[a,b]]));
}
for(a=0,7, for(b=0,7, ajoute_si_primitif(a,b)));
print("nb paires primitives mod 8 : ", #prims);

vues=[]; nclasses=0;
ajoute_classe(p)={
  my(c=classe_de(p));
  my(deja=0);
  for(j=1,#vues, if(vues[j]==c, deja=1));
  if(!deja, vues=concat(vues,[c]); nclasses=nclasses+1);
}
for(i=1,#prims, ajoute_classe(prims[i]));
print("nb classes distinctes : ", nclasses);
affiche_rep(j)={ print("  ", vues[j][1]); }
for(j=1,#vues, affiche_rep(j));

mes12=[[0,1],[1,1],[2,1],[3,1],[4,1],[5,1],[6,1],[7,1],[1,0],[1,4],[1,2],[3,2]];
print("\nmes 12 points, classes:");
mesclasses=[];
calcule_mesclasse(i)={ mesclasses=concat(mesclasses,[classe_de(mes12[i])]); }
for(i=1,#mes12, calcule_mesclasse(i));

ndoublons=0;
teste_doublon(i,j)={
  if(mesclasses[i]==mesclasses[j], ndoublons=ndoublons+1; printf("  DOUBLON entre (%d,%d) et (%d,%d)\n", mes12[i][1],mes12[i][2],mes12[j][1],mes12[j][2]));
}
for(i=1,#mesclasses, for(j=i+1,#mesclasses, teste_doublon(i,j)));
print("doublons trouves : ", ndoublons);

ncouvert=0;
teste_couverture(j)={
  my(trouve=0);
  for(i=1,#mesclasses, if(mesclasses[i]==vues[j], trouve=1));
  if(trouve, ncouvert=ncouvert+1);
}
for(j=1,#vues, teste_couverture(j));
print("classes couvertes par mes 12 points : ", ncouvert, " / ", nclasses);
quit
