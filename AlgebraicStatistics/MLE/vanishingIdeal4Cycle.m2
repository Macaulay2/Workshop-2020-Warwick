restart
loadPackage "GraphicalModels"
G=graph{{1,2},{2,3},{3,4},{1,4}}
R=gaussianRing G
I=gaussianVanishingIdeal(R)
betti trim I

flatten entries gens gb I

 netList (trim I)_*
