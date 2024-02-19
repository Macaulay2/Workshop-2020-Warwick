R=QQ[l_1..l_6,m_1..m_4,a_1..a_4]
M=matrix{{l_1,l_4,0,l_4},{l_4,l_2,l_5,0},{0,l_5,l_3,l_6},{l_4,0,l_6,l_2}}
N=matrix{{0,m_2,m_3,-m_2},{m_2,m_1,0,m_4},{m_3,0,0,0},{-m_2,m_4,0,-m_1}}
X=matrix {{a_1},{a_2},{a_3},{a_4}}
M*X

X*transpose(X)+N

leadingPrincipalMinors=M->(
mm:=numcols M;
ss:=toList apply(0..(numcols M-1),i->toList(0..i));
for s in ss list M_s^s
);

L=leadingPrincipalMinors(X*transpose(X)+N)
netList apply(L,i->factor det i)

restart
R=QQ[s_(1,1)..s_(1,4),s_(2,2)..s_(2,4),s_(3,3),s_(3,4),s_(4,4)]
S=genericSymmetricMatrix(R,s_(1,1),4)
I2=minors(3,S);
dim I2,codim I2, degree I2

X=random(QQ^4,QQ^2);              
M=X*transpose(X);
l1=X_{0}
l2=X_{1}
X_{0}*(transpose(X_{0}))+X_{1}*(transpose(X_{1}))==M
rank jacobian(I2)
I2s=sub(jacobian(I2),matrix {flatten toList apply(0..3,i->toList apply(i..3,j->M_(i,j)))})
rank I2s
Ps=matrix {apply(gens R,flatten toList apply(0..3,i->toList apply(i..3,j->M_(i,j))),(i,j)->i-j)}

J2=trim ideal flatten entries (Ps*I2s);
betti (trim J2)
dim J2, codim J2, degree J2
codim J2==rank I2s


I1=minors(2,S);
dim I1,codim I1, degree I1
X1=random(QQ^4,QQ^1);              
l=X1
M1=X1*transpose(X1);
rank jacobian(I1)
jacobian I1
I1s=sub(jacobian(I1),matrix {flatten toList apply(0..3,i->toList apply(i..3,j->M1_(i,j)))})
rank I1s
Ps1=matrix {apply(gens R,flatten toList apply(0..3,i->toList apply(i..3,j->M1_(i,j))),(i,j)->i-j)}

J1=trim ideal flatten entries (Ps1*I1s);
betti J1
dim J1, codim J1, degree J1
codim J1==rank I1s

--can we say that T_M V^4_1 is the degree 2 part of the ideal spanned by l in R?
