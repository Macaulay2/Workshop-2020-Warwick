restart
loadPackage("GraphicalModelsMultiTrek",Reload=>true)
--G = mixedGraph(digraph {{b,{c,d}},{c,{d}}},bigraph {{a,d}})
g = mixedGraph(digraph {{b,d},{c,d},{b,c}},bigraph {{a,d}})
R = gaussianRing G
T = trekSeparation G
T = apply(T,s -> {sort s#0, sort s#1, sort s#2, sort s#3})
L = {{{a}, {b, c}, {}, {}}, {{b, c}, {a, b}, {}, {b}}, {{a, b}, {b, c}, {}, {b}}, {{b, c}, {a, c}, {}, {c}}, {{b, c}, {a, d}, {}, {d}}}
sort T=== sort L

multiTrekSeparation(g,4)
((set subsets vertices G)^**3)/splice
toList(((set subsets vertices G)^**3)/splice)

toList(((set subsets vertices G)^**3)/splice)

G=graph collateVertices g;
(graph G#Digraph)#b

restart
loadPackage("GraphicalModelsMultiTrek",Reload=>true)
g=digraph {{b,d},{c,d},{b,c}}
tops(g,2,{{b},{c,d}},{{c},{c}})


--
k=2
Slist = {{b},{c,d}}
Alist = {{c},{b,c}}
topList := {};
    Glist:=apply(Alist,A->deleteVertices(g,A));
    vert := sort vertices g;
    for v in vert do(
	isTop := true;
	i:=0;
	descd := {};
	while (isTop and i<k) do(
	    pathExists := false;
	    descd = toList(descendants(Glist#i,v));
	    print descd;
	    ddd := 0;
	    print "a";
	    while ((not pathExists) and ddd<length(descd)) do(
		if member(descd#ddd,Slist#i) then(pathExists=true;);
		ddd=ddd+1;
		);
	    print "b";
	    if not pathExists then (isTop=false;);
	    i=i+1; 
	    print "c";
	    );
	print "d";
	if isTop then(topList=append(topList,v););
	print "hi";	