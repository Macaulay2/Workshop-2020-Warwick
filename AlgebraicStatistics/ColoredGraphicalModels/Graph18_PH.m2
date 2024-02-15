
restart
L1=value get "M2_rk1Data.txt";
apply(L1,i->rank(i))
apply(0..9,i->eigenvalues L1_i)

L2=value get "M2_rk2Data.txt";
apply(L2,i->rank(i))
apply(0..9,i->eigenvalues L2_i)


L3=value get "M2_rk3Data.txt";
apply(L3,i->rank(i))
apply(0..9,i->eigenvalues L3_i)
eigenvalues L1_0
