-- implement syzLLL

-- Strategy I: Degree by degree
-- Strategy II: Use lower degrees to simplify higher degrees 
-- (a,b)


restart
--installPackage"RandomObjects"
uninstallPackage"RandomObjects"
installPackage"RandomObjects"

uninstallPackage"RandomSpaceCurves"
installPackage"RandomSpaceCurves"



help loadPackage

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
    -- compute nunber of different degrees
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

-- apply LLL to a matrix with several degrees
matrixLLL = method()

matrixLLL (Matrix) := (M) -> (
    degs := unique degrees source M;
    MLLL := fold((a,b)->(a|b),
        apply(degs,d->matrixLLLoneDegree selectDegree(M,d))
        );
    -- copy degrees from M to MLLL
    -- 
    -- CAVEAT: this will not work if the degrees in M 
    --         are mixed up
    MLLLhomogeneous := 0*M+MLLL;
    -- TEST: is it really honmogeneous
    assert isHomogeneous MLLLhomogeneous;
    return MLLLhomogeneous
    )

-- apply LLL to a syzygy matrix Degree by Degree

syzLLL1 = method()
syzLLL1 (Matrix) := M -> matrixLLL syz M

-- TEST: does this indeed give syzygies?
assert (
    M := random(RZZ^{3:0},RZZ^{9:-1});
    0 == M * syzLLL1 M
    )

-- maximum coefficient of a matrix
height (Matrix) := M -> lift(max flatten entries last coefficients M,ZZ)

-- one example
tally apply(10,i->(
    h := 3;
    print "";
    time betti (M = random(RZZ^{3:0},RZZ^{9:-1},Height=>h));
    time betti (syzM = syzLLL1 M);
    --degs := unique degrees source syzM;
    --apply(degs,d->height matrixLLLoneDegree selectDegree(syzM,d))
    apply(rank source syzM,i->height (syzM_{i}))  
    ))
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

-- Strategy II

betti syzM
super 

basis(1,ker M)
mons1 = super basis(1,RZZ)


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
height gens I

time betti 
time codim singIQQ

