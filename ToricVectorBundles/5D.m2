 --Bu dosyada P(1,w1,...,wn) tipindeki WPS'ler calisiliyor
 --normalToricvarieties paketinde WPS ler de var...
 restart
 F=ZZ/3;
 S=F[x_0..x_5, Degrees => {1,1,1,1,1,2}];
 J=ideal(x_0^2-x_1^2,x_0^2-x_2^2,x_0^2-x_3^2,x_0^2-x_4^2,x_0^4-x_5^2);
 J=ideal(x_0^2-x_1^2,x_0^2-x_2^2,x_0^2-x_3^2,x_0^2-x_4^2,x_0^2-x_5);
 B=ideal(x_0,x_1,x_2,x_3,x_4,x_5);
 IY=saturate (J,B);
IY==J
----Asagida HFsi 32 den farkli olan k lar elde ediliyor. Elde edilen son k=a_Y
k=0;
while hilbertFunction(k+1,J) != 32
do k=k+1;

--dimensions are
for i to k-1
 list hilbertFunction(i+1,J)
----
for i to k-1
 list transpose basis({i+1},S/J);--monomların matrisi
-- 
 m=oo;
-----bazı komutlar
 oo / print @@ print;--lists the previous output
for i to 4
 list m_0_(i,0)-----Matrisin elemanlarindan liste yapmak
Matrix ^ List -- select rows
Matrix ^ ZZ -- power
Matrix _ List -- select columns
submatrix(f, {1,2}, {3,4})--f nin 1 ve 2. satiri ile 3 ve 4. sutunu alinir
scan(BasicList,Function) -- apply a function to each element of a list
positions(VisibleList,Function) -- which elements of a list satisfy a condition
-----------T nin elemanlarini kartezyen carpim ile bulsak kume veriyor liste degil
S=set{1,2};
s=set{1};
T=(s ** S^**5)/splice ;
Ts=subsets(T,31); --T nin 31 elemanli altkumelerinin kumesi
--
for i to #Ts-1
list toList Ts_i; --T kumesini listeye cevirir
K=oo;
T=apply(#Ts-1,i-> apply (#K, j-> toList K_i_j));

-----------------Bu komut daha iyi
loadPackage "RationalPoints";
 R = ZZ/3[x,y,z,t,u];
 P=rationalPoints ideal(x^(2)-1,y^2-1,z^2-1,t^(2)-1,u^2-1); --points in 5-torus 
T= apply(P,p-> flatten {1,p});--points in 5-torus in homogeneous coordinates
Ts=subsets(T,31);
--------Generating matrix of the toric code for degree (1,2)
SentGap:=(T,m)->(
    a=apply(#T,j->  sub( m, matrix {T_j}));
    for j to (numgens target m)-1 
    do (
	if j==0 then c=[]; 
	b=[]; 
	for i to #a-1 
	do (
	    b=append(b, (a_i_(j,0))); 
	    ) ;  
	c=append(c,b); 
	); 
   f="Dropbox/m2-dosyalar/deneme"<< "M:="<< c <<";" << endl<<close ;
    ) ;

SentGap(T,m_0)
value get "deneme" --Bu dosya GAP ta matris tanimlamaya uygun, Macaulayda array.
--GAP
Downloads/gap4r7/bin/gap.sh
LoadPackage( "guava" );
Read("Dropbox/m2-dosyalar/deneme"); (Bu macaulaydan gelen dosyayi okutuyor GAPa)
C := GeneratorMatCode(M,GF(3));
MinimumDistance(C);

 