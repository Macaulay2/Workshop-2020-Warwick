-- implement syzLLL

-- Strategy I: Degree by degree
-- Strategy II: Use lower degrees to simplify higher degrees 
-- (a,b)


--restart
--installPackage"RandomObjects"
--uninstallPackage"RandomObjects"
--installPackage"RandomObjects"

--uninstallPackage"RandomSpaceCurves"
--installPackage"RandomSpaceCurves"

--help loadPackage

restart
RQQ = QQ[x_0..x_3]
RZZ = ZZ[x_0..x_3]


-- select columns of a given degree
--
-- Input:
--  M: a matrix
--  d: a degree
selectDegree = method()

selectDegree (Matrix,List) := (M,d) -> (
    degs := degrees source M;
    columnsOK := select(rank source M,i->d == degs#i);
    M_columnsOK
    )

-- apply LLL to columns of one degree of a Matrix
matrixLLLoneDegree = method()

matrixLLLoneDegree (Matrix) := M -> (
    -- compute number of different degrees
    numDegs := #unique degrees source M;
    -- TEST: does the source of M have only one degree?
    if numDegs != 1 then error "Matrix has multiple degrees";
    -- Do LLL on the coefficients
    (mons,coeffs) := coefficients M;
    LLLM = mons * LLL lift(coeffs,ZZ);
    -- TEST: does this still define the same cokernel?
    assert (coker M == coker LLLM);
    return LLLM
    )

-- TEST: does this work with syzygies?
assert (
    M := random(RZZ^{3:0},RZZ^{9:-1});
    syzM := syz M;
    0 == M * matrixLLLoneDegree selectDegree(syzM,{2})
    )


matrixLLL = method(Options=>{strategyI=>false})

matrixLLL (Matrix) := Matrix => opt -> (M) -> (
    R := ring M; 
    degs := unique degrees source M;
    if opt.strategyI then (
    -- Stategy I: apply LLL to a matrix with several degrees	
    MLLLinhomogeneous := fold((a,b)->(a|b),
        apply(degs,d->matrixLLLoneDegree selectDegree(M,d))
        );
    -- copy degrees from M to MLLL
    -- CAVEAT: this will not work if the degrees in M 
    --         are mixed up
    MLLL := 0*M+MLLLinhomogeneous;)
    -- Stategy II: apply LLL to basis of columns in several degrees
    else ( 
        MLLL = fold((a,b)->(a|b),
        apply(degs,d-> (Md := matrixLLLoneDegree (M * matrix basis(d,source M));
	    map(target M, R^{rank source Md: -d}, Md)))
        ););
    -- TEST: is it really honmogeneous
    assert isHomogeneous MLLL;
    return MLLL
    )


matrixLLL2 = method()

matrixLLL2 (Matrix) := (M) -> (
    R := ring M;
    degs := unique degrees source M;
    MLLL := fold((a,b)->(a|b),
        apply(degs,d-> (Md := matrixLLLoneDegree (M * matrix basis(d,source M));
	    map(target M, R^{rank source Md: -d}, Md)))
        );
    -- TEST: is it really honmogeneous
    assert isHomogeneous MLLL;
    return MLLL
    )

-- apply LLL to a syzygy matrix Degree by Degree

syzLLL1 = method(Options=>{strategyI=>false})
syzLLL1 (Matrix) := Matrix => opt -> M -> matrixLLL(syz M, strategyI=>opt.strategyI)

syzLLL2 = method()
syzLLL2 (Matrix) := M -> matrixLLL2 syz M

-- TEST: does this indeed give syzygies?
assert (
    M = random(RZZ^{3:0},RZZ^{9:-1});
    0 == M * syzLLL1 M
    )

assert (
    M = random(RZZ^{3:0},RZZ^{9:-1});
    0 == M * syzLLL2 M
    )

-- maximum coefficient of a matrix
height (Matrix) := M -> lift(max flatten entries last coefficients M,ZZ)

end;

-- experiment 1: comparing Stategy I and II
restart;
load "syzLLLtests.m2"

h = 3
apply(10, i-> (
	betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h));
	print("StrategyI:  "| toString (time height syzLLL1(M, strategyI=>true)));
	print("StrategyII: "| time height syzLLL1 M);
	print"";	
	));

betti res coker M
betti res coker sub(M,RQQ)
syzM = syz sub(M,RQQ);
N = syzM * random( source syzM, RQQ^{6:-2,1:-3});
I = ideal syz N;
dim I, degree I, genus I
betti (resI = res I)
ann coker transpose resI.dd_3 == ann coker sub(M,RQQ)

time tally apply(20, i->(
    h := 10;
    betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h));
    syzM := syzLLL1 M;
    (hsyz := height syzM, log(hsyz)/log(h))
    ))

-*
height = 2
-- used 2.43878 seconds
      Tally{(3, 1.58496) => 1    }
            (4, 2) => 2
            (5, 2.32193) => 3
            (6, 2.58496) => 5
            (7, 2.80735) => 1
            (13, 3.70044) => 1
            (19, 4.24793) => 1
            (20, 4.32193) => 2
            (22, 4.45943) => 1
            (58, 5.85798) => 1
            (176, 7.45943) => 1
            (11116, 13.4404) => 1

height = 5
     -- used 7.64966 seconds
      Tally{(456, 3.80412) => 1 }                                                              
            (561, 3.93288) => 1
            (616, 3.99099) => 1
            (673, 4.04597) => 1
            (677, 4.04966) => 2
            (789, 4.14478) => 1
            (813, 4.1634) => 1
            (850, 4.19105) => 1
            (855, 4.1947) => 1
            (867, 4.20336) => 1
            (986, 4.28327) => 1
            (1148, 4.37779) => 1
            (1468, 4.53056) => 1
            (1628, 4.59484) => 1
            (2666, 4.9013) => 1
            (2794, 4.93044) => 1
            (485054689, 12.4266) => 1
            (23618346658494, 19.1328) => 1
            (2689533985020386..., 409.788) => 1	  
	    
height = 10
-- used 30.1941 seconds  
o21 = Tally{(14096, 4.1491) => 1    }                                                         
            (23048, 4.36263) => 1
            (23507, 4.3712) => 1
            (23541, 4.37182) => 1
            (27502, 4.43936) => 1
            (28769, 4.45892) => 1
            (32097, 4.50646) => 1
            (35919, 4.55532) => 1
            (41462, 4.61765) => 1
            (41891, 4.62212) => 1
            (42149, 4.62479) => 1
            (44614, 4.64947) => 1
            (45504, 4.65805) => 1
            (55490, 4.74421) => 1
            (286954, 5.45781) => 1
            (356186, 5.55168) => 1
            (364816, 5.56207) => 1
            (380215, 5.58003) => 1
            (756432, 5.87877) => 1
            (28824273434..., 129.46) => 1
*-	    

time tally apply(20, i->(
    h := 10;
    betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h));
    syzM := syzLLL2 M;
    (hsyz := height syzM, log(hsyz)/log(h))
    ))

-*
height = 2
-- used 1.50505 seconds
      Tally{(2, 1) => 1       }
            (3, 1.58496) => 4
            (4, 2) => 5
            (5, 2.32193) => 3
            (6, 2.58496) => 3
            (7, 2.80735) => 2
            (8, 3) => 1
            (10, 3.32193) => 1
	    
height = 5
-- used 10.1318 seconds
      Tally{(414, 3.74408) => 1 }
            (549, 3.91944) => 1
            (655, 4.02913) => 1
            (667, 4.04041) => 1
            (675, 4.04782) => 1
            (761, 4.12233) => 1
            (783, 4.14004) => 1
            (802, 4.15493) => 1
            (862, 4.19976) => 1
            (963, 4.2686) => 1
            (974, 4.27566) => 1
            (983, 4.28138) => 1
            (1041, 4.317) => 1
            (1197, 4.40376) => 1
            (1221, 4.41609) => 1
            (1230, 4.42065) => 1
            (1297, 4.45361) => 2
            (1391, 4.49708) => 1
            (1474, 4.5331) => 1	  

height = 10
-- used 25.6716 seconds
      Tally{(15352, 4.18616) => 1}
            (21221, 4.32677) => 1
            (24900, 4.3962) => 1
            (25273, 4.40266) => 1
            (31110, 4.4929) => 1
            (31932, 4.50423) => 1
            (32216, 4.50807) => 1
            (36091, 4.5574) => 1
            (38358, 4.58386) => 1
            (38534, 4.58584) => 1
            (38917, 4.59014) => 1
            (42095, 4.62423) => 1
            (43328, 4.63677) => 1
            (43549, 4.63898) => 1
            (44486, 4.64822) => 1
            (44962, 4.65285) => 1
            (55256, 4.74238) => 1
            (61134, 4.78628) => 1
            (76772, 4.8852) => 1
            (82249, 4.91513) => 1	      
*-	     
 

-- one example
-*
tally apply(10,i->(
    h := 4;
    print "";
    time betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h));
    time betti (syzM = syzLLL1 M);
    --degs := unique degrees source syzM;
    --apply(degs,d->height matrixLLLoneDegree selectDegree(syzM,d))
    apply(rank source syzM,i->height (syzM_{i}))  
    ))
*-
-- the height of the result seems to be about h^4
-- but some abnormally big coefficients appear (
-- 
-- 549 => 1            
-- 609 => 1
-- 706 => 1
-- 857 => 1
-- 955 => 1
-- 1914 => 1
-- 1972 => 1
-- 11112 => 1
-- 67471871 => 1
-- 121778925236300 => 1

-- these might come from non optimally chosen
-- bases for higher degree syzygies.


restart;
load "syzLLLtests.m2"
loadPackage "FastLinAlg"
loadPackage("PruneComplex",Reload=>true)

RZZ = ZZ[x_0..x_3]
RQQ = QQ[x_0..x_3]

h = 2
betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h));
betti (syzM = syz M);
LLLdeg2 = matrixLLL(syzM * basis(2,source syzM));
LLLdeg3 = matrixLLL(syzM * basis(3,source syzM));
betti (N = LLLdeg2|LLLdeg3)
betti (syzN = syz N)

-- we want to minimize LLLdeg3

-- using pruneComplex
C = chainComplex({N,syzN})
betti pruneDiff(C,1)
C' = pruneComplex (C)
betti (C'.dd_2)
max flatten entries last coefficients (N)
max flatten entries last coefficients (C'.dd_1)
assert (M * C'.dd_1 == 0)
-- we get a smaller matrix 
-- and the height of the coefficients does not change

-- using FastLinAlg -> getSubmatrixOfRank
betti (constantSyz = sub(syzN_{0..31}^{6..43},RQQ))
rkConstantSyz = rank constantSyz;
subList = getSubmatrixOfRank(rkConstantSyz, constantSyz) 
rowsDeleted = apply(subList#0, i-> i + 6)
betti (improvedN = submatrix'(N,,rowsDeleted))
N' = improvedN * random(source improvedN, RZZ^{6:-2,1:-3},Height=>2);
betti (I = ideal matrixLLL syz N')
IQQ = sub(I,RQQ);
singIQQ = ideal singularLocus IQQ;
hI = height gens I;
result =  (
    hI,
    log(hI)/log(h),
    codim singIQQ
    )




-- experiment without minimizing: 

h = 2
time tally apply(20,i->(
        betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h));
        betti (syzM = syz M);
        LLLdeg2 = matrixLLL(syzM * basis(2,source syzM));
        LLLdeg3 = matrixLLL(syzM * basis(3,source syzM));
        -- strategy II is much better.
        N = (
            LLLdeg2*random(source LLLdeg2,RZZ^{6:-2},Height=>h)|
            LLLdeg3*random(source LLLdeg3,RZZ^{1:-3},Height=>h)
            );
        betti (syzN = syz N);
        I = ideal matrixLLL gens ideal syzN;
        IQQ = sub(I,RQQ);
        if codim IQQ == 2 then (
            singIQQ = ideal singularLocus IQQ;
            hI = height gens I;
            result =  (
                hI,
                log(hI)/log(h),
                --betti res IQQ,
                codim singIQQ
                )
            ) else result = {};
        return result
        )
    )

 -- used 9.11849 seconds (for h=2)

-- (34, 5.08746, 3) => 1  }
-- (235, 7.87652, 3) => 1
-- (368, 8.52356, 4) => 1
-- (393, 8.61839, 4) => 1
-- (552, 9.10852, 4) => 1
-- (962, 9.90989, 4) => 1
-- (1056, 10.0444, 3) => 1
-- (1363, 10.4126, 4) => 1
-- (2454, 11.2609, 4) => 1
-- (4017, 11.9719, 4) => 1
-- {} => 10

h = 2
time tally apply(20,i->(
        M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h);
        betti (syzM = syz M);
        betti (N = syzM * random(source syzM, RZZ^{6:-2,1:-3},Height=>h));
        betti (I = ideal syz N);
        IQQ = sub(I,RQQ);
        if codim IQQ == 2 then (
            singIQQ = ideal singularLocus IQQ;
            hI = height gens I;
            result =  (
                hI,
                log(hI)/log(h),
                --betti res IQQ,
                codim singIQQ
                )
            ) else result = {};
        return result
    ))
-- used 4.92308 seconds

-- for h=2
--
-- (635, 9.31061, 4) => 1   }
-- (1142, 10.1573, 4) => 1
-- (1757, 10.7789, 4) => 1
-- (2826, 11.4645, 4) => 1
-- (2921, 11.5122, 3) => 1
-- (20572, 14.3284, 4) => 1
-- (21519, 14.3933, 3) => 1
-- (60927, 15.8948, 4) => 1
-- (72317, 16.142, 4) => 1
-- (828586, 19.6603, 4) => 1
-- {} => 10


restart;
RZZ = ZZ[x_0..x_3]
h = 2
betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h))
betti (syzM = syz M)
betti (M1 = syzM * basis(2,source syzM))
betti (M2 = syzM * basis(3,source syzM))
betti (syz (M1|M2))

methods random
code 10

