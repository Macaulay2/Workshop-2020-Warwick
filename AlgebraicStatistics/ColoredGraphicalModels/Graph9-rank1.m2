-- GRAPH 9
restart
loadPackage "EigenSolver"
n=8
R=QQ[l_1..l_n,t_1..t_n,s_11..s_14,s_22..s_24,s_33..s_34,s_44]
R=QQ[s_11..s_14,s_22..s_24,s_33..s_34,s_44,t_8,t_7,t_6,t_5,t_4,t_3,t_2,t_1]
--K=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
S=matrix{{s_11,s_12,s_13,s_14},{s_12,s_22,s_23,s_24},{s_13,s_23,s_33,s_34},{s_14,s_24,s_34,s_44}}--estimate

--sufficient statistics
varList=flatten {toList(s_11..s_14), toList(s_22..s_24),s_33,s_34,s_44}
Istat=ideal(t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34,t_7-2*s_13,t_8-2*s_24)
Istatorig=ideal(t_1-s_11,t_2-s_22-s_44,t_3-s_33,t_4-2*(s_12+s_14),t_5-2*s_23,t_6-2*s_34)

I_G1=eliminate(varList,minors(2,S)+Istat)
I_G1orig=eliminate(varList,minors(2,S)+Istatorig)
J=minors(2,S)+Istat

eliminate({t_7,t_8},I_G1)==I_G1orig

statCand={738245279/15389164,
      540127267/10939086,
      363531347/20054991,
      2110596472/21603499,
      1167167116/27015947,
      568933285/33760698}
 
P1={t_1 => 55, t_2 => 329/4, t_3 => 35, t_4 => 365/2, t_5 => 68, t_6 => 109/2}
for i from 0 to 1 do print sub(J_i,P1)
v1=vector{55,329/4,35,365/2,68,109/2}
P2={t_1 => 25, t_2 => 277/4, t_3 => 131/4, t_4 => 95, t_5 => 111/2, t_6 => 58}
for i from 0 to 1 do print sub(J_i,P2)
v2=vector{25,277/4,131/4,95,111/2,58}
-- Point in the intersection of the segment between P1 and P2 and relevant
test={44.4548, 77.6804, 34.2091, 151.743, 63.6062, 55.7303}  

statCand={44.4548, 77.6804, 34.2091, 151.743, 63.6062, 55.7303}
statCand=apply(statCand,i->lift(i,QQ))


dim I_G1  
  
partialSol=sub(I_G1,{t_1 => statCand_0, t_2 => statCand_1, t_3 => statCand_2, t_4 => statCand_3, t_5 => statCand_4, t_6 => statCand_5})  
sub((gens partialSol)_1,RR)

dim partialSol

entries (gens partialSol)_1
entries (gens partialSol)_1

partialSol6=ideal (flatten entries (gens partialSol)_6)

sols6=zeroDimSolve sub(partialSol6,QQ[t_7])

partialSol0=ideal (flatten entries (gens partialSol)_0)
sols0=zeroDimSolve sub(partialSol6,QQ[t_8])

t_8=lift(.365332,QQ)

partialSol7=sub(partialSol,{t_8 => lift(.365332,QQ)})
partialSol7=sub(partialSol,{(support partialSol)_1 => lift(.365332,QQ)})
sols7=zeroDimSolve sub(partialSol7,QQ[t_7])
--in the model??

partialSol= ideal(-(727062694/20054991)*t
      _8+332020110724928030/
      456038613925503,-2*t_7*t_8-(
      540127267/5469543)*t_7+
      57803589404431331843249668/
      9852029739900990134997,-(
      2110596472/21603499)*t_7+
      40437125289647364115333777/
      7018053020032249449492,-(
      54774647992902863/
      912077227851006)*t_7+
      1534535956879215568/
      433257978013509,
      14494848254481264849641457689
      358703139822487799/
      10138979111771840723712151182
      649964114440452,-t_7^2+
      268375300691260813/
      77157386379381,(738245279/
      3847291)*t_8-
      13769489097813917723384588388
      91/19641934907685506606230234
      026)
  
 a=(332020110724928030/456038613925503)/(727062694/20054991)
 
 a'=(1376948909781391772338458838891/19641934907685506606230234026)/(738245279/3847291)

sub(a'-a,RR)
