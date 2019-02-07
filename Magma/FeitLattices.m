K:=QuadraticField(-3);
R:=Integers(K);
KK:=FieldOfFractions(R);

TrueBasis:=function(L)
 PB:=PseudoBasis(L);
 res:=[];
 for v in PB do
  _,g:=IsPrincipal(v[1]);
  Append(~res,g*v[2]);
 end for;
 return res;
end function;

AnLat:=function(n)
 mat:=RMatrixSpace(R,n+1,1)![1: i in [1..n+1]];
 B:=ChangeRing(BasisMatrix((Kernel(mat))),K);
 return B;
end function;

DnLat:=function(n,alpha)
 mat:=RMatrixSpace(R,n,1)![1: i in [1..n]];
 B:=BasisMatrix(Kernel(mat));
 B:=VerticalJoin(B,MatrixRing(R,n)!alpha);
 B:=Matrix(TrueBasis(Module(Rows(B))));
 B:=ChangeRing(B,K);
 return B;
end function;

omega:=R!((K.1-1)/2);


B:=VerticalJoin(AnLat(5),1/K.1*KMatrixSpace(K,1,6)![1,omega,omega^2,1,omega,omega^2]);
U5bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
U5:=HermitianLattice(U5bas*HermitianTranspose(U5bas));

B:=VerticalJoin(DnLat(6,K.1),1/K.1*KMatrixSpace(K,1,6)![1,1,1,1,1,1]);
U6bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
U6:=HermitianLattice(U6bas*HermitianTranspose(U6bas));


//In Lati we will store the indecomposable unimodular lattices in dimension i according to table I in Feit

Lat6:=[U6];  

E8:=HermitianLattice(ChangeRing(CartanMatrix(RootSystem("E8")),K));
Lat8:=[E8];

B:=VerticalJoin(DnLat(9,K.1),1/K.1*KMatrixSpace(K,1,9)![1,1,1,1,1,1,1,1,1]);
lat9bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat9:=HermitianLattice(lat9bas*HermitianTranspose(lat9bas));
Lat9:=[lat9];

B:=VerticalJoin(DnLat(10,2),1/2*KMatrixSpace(K,1,10)![1+2*omega,1,1,1,1,1,1,1,1,1]);
lat101bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat101:=HermitianLattice(lat101bas*HermitianTranspose(lat101bas));

U5plusU5:=VerticalJoin(HorizontalJoin(U5bas,KMatrixSpace(K,5,6)!0),HorizontalJoin(KMatrixSpace(K,5,6)!0,U5bas));
B:=VerticalJoin(U5plusU5,1/(2*K.1)*KMatrixSpace(K,1,12)![1,1,1,1,1,-5,1,1,1,1,1,-5]);
lat102bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat102:=HermitianLattice(lat102bas*HermitianTranspose(lat102bas));
Lat10:=[lat101,lat102];

B:=VerticalJoin(AnLat(11),1/(2*K.1)*KMatrixSpace(K,1,12)![1,1,1,1,1,1,1,1,1,1,1,-11]);
lat111bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat111:=HermitianLattice(lat111bas*HermitianTranspose(lat111bas));

E6:=HermitianLattice(ChangeRing(CartanMatrix(RootSystem("E6")),K));
d5:=DnLat(5,K.1);
d5lat:=HermitianLattice(d5,MatrixRing(K,5)!1);
DE6:=Dual(E6);
alpha:=[1,0,-1,0,1,-1];
Z:=DirectSum(E6,d5lat);
lat112bas:=ChangeRing(Matrix(TrueBasis(Module(Generators(Module(Z)) cat [ChangeRing(Vector([1/3*x: x in alpha] cat [1/K.1: i in [1..5]]),KK)]))),K);
lat112:=HermitianLattice(lat112bas*Z`Form*HermitianTranspose(lat112bas));
Lat11:=[lat111,lat112];

II:=func<n| HermitianLattice(MatrixRing(K,n)!1)>;
Unimod11:=[II(11),DirectSum(Lat6[1],II(5)),DirectSum(Lat8[1],II(3)),DirectSum(Lat9[1],II(2)),DirectSum(Lat10[1],II(1)),DirectSum(Lat10[2],II(1)),Lat11[1],Lat11[2]];


//Now for the lattices in dimension 12, according to table II in Feit
A6plusA6:=VerticalJoin(HorizontalJoin(AnLat(6),KMatrixSpace(K,6,7)!0),HorizontalJoin(KMatrixSpace(K,6,7)!0,AnLat(6)));
B:=VerticalJoin(A6plusA6,1/7*KMatrixSpace(K,1,14)!([1,1,1,1,1,1,-6] cat [3*K.1*x: x in [1,1,1,1,1,1,-6]]));
lat121bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat121:=HermitianLattice(lat121bas*HermitianTranspose(lat121bas));

B:=VerticalJoin(DnLat(12,K.1),1/K.1*KMatrixSpace(K,1,12)![1: i in [1..12]]);
lat122bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat122:=HermitianLattice(lat122bas*HermitianTranspose(lat122bas));

B:=VerticalJoin(DnLat(12,2),1/2*KMatrixSpace(K,1,12)![1: i in [1..12]]);
lat123bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat123:=HermitianLattice(lat123bas*HermitianTranspose(lat123bas));

B:=VerticalJoin(AnLat(12),1/(4*omega+1)*KMatrixSpace(K,1,13)!([1: i in [1..12]] cat [-12]));
lat124bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat124:=HermitianLattice(lat124bas*HermitianTranspose(lat124bas));

lat125bas:=Parent(lat124bas)![Conjugate(x): x in Eltseq(lat124bas)];
lat125:=HermitianLattice(lat125bas*HermitianTranspose(lat125bas));

E7:=HermitianLattice(ChangeRing(CartanMatrix(RootSystem("E7")),K));
Z:=DirectSum(E7,U5);
alpha:=Vector([K!0,1,0,0,1,0,1]);
xx:=Solution(U5bas,1/(2*K.1)*Vector([K!1,1,1,1,1,-5]));
lat126bas:=ChangeRing(Matrix(TrueBasis(Module(Generators(Module(Z)) cat [ChangeRing(Vector(Eltseq(1/2*alpha) cat Eltseq(xx)),KK)]))),K);
lat126:=HermitianLattice(lat126bas*Z`Form*HermitianTranspose(lat126bas));


B:=VerticalJoin(HorizontalJoin(AnLat(8),KMatrixSpace(K,8,4)!0),HorizontalJoin(KMatrixSpace(K,4,9)!0,DnLat(4,K.1)));
B:=VerticalJoin(B,KMatrixSpace(K,1,13)!([1/(3*K.1)*x: x in [1,1,1,1,1,1,1,1,-8]] cat [1/K.1*x: x in [1,1,1,1]]));
lat127bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat127:=HermitianLattice(lat127bas*HermitianTranspose(lat127bas));

B:=VerticalJoin(HorizontalJoin(DnLat(6,2),KMatrixSpace(K,6,6)!0),HorizontalJoin(KMatrixSpace(K,6,6)!0,DnLat(6,2)));
B:=VerticalJoin(B,KMatrixSpace(K,1,12)!([1/2*x: x in [1,1,1,1,1,2*omega+1]] cat [0,0,0,0,0,1]));
B:=VerticalJoin(B,KMatrixSpace(K,1,12)!([0,0,0,0,0,1] cat [1/2*x: x in [1,1,1,1,1,2*omega^2+1]]));
lat128bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat128:=HermitianLattice(lat128bas*HermitianTranspose(lat128bas));

B:=VerticalJoin(HorizontalJoin(DnLat(6,K.1),KMatrixSpace(K,6,6)!0),HorizontalJoin(KMatrixSpace(K,6,6)!0,DnLat(6,K.1)));
B:=VerticalJoin(B,KMatrixSpace(K,1,12)!([1/K.1*x: x in [1,1,1,1,1,1]] cat [0,0,0,0,0,omega^2]));
B:=VerticalJoin(B,KMatrixSpace(K,1,12)!([0,0,0,0,0,omega] cat [1/K.1*x: x in [1,1,1,1,1,1]]));
lat129bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat129:=HermitianLattice(lat129bas*HermitianTranspose(lat129bas));

B:=HorizontalJoin(DnLat(3,K.1),KMatrixSpace(K,3,9)!0);
B:=VerticalJoin(B,HorizontalJoin(KMatrixSpace(K,3,3)!0,HorizontalJoin(DnLat(3,K.1),KMatrixSpace(K,3,6)!0)));
B:=VerticalJoin(B,HorizontalJoin(KMatrixSpace(K,3,6)!0,HorizontalJoin(DnLat(3,K.1),KMatrixSpace(K,3,3)!0)));
B:=VerticalJoin(B,HorizontalJoin(KMatrixSpace(K,3,9)!0,DnLat(3,K.1)));
d3:=DnLat(3,K.1);
d3lat:=HermitianLattice(d3*HermitianTranspose(d3));
d3d:=Dual(d3lat);
gen:=Generators(quo<Module(d3d)|Module(d3lat)>);
for xx in gen do
 for yy in gen do
  x:=xx*ChangeRing(d3,KK);
  y:=yy*ChangeRing(d3,KK);
  B:=VerticalJoin(B,KMatrixSpace(K,1,12)![Eltseq(x) cat Eltseq(y) cat Eltseq(x+y) cat Eltseq(x-y)]);
 end for;
end for; 
lat1210bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat1210:=HermitianLattice(lat1210bas*HermitianTranspose(lat1210bas));

B:=HorizontalJoin(DnLat(4,2),KMatrixSpace(K,4,8)!0);
B:=VerticalJoin(B,HorizontalJoin(KMatrixSpace(K,4,4)!0,HorizontalJoin(DnLat(4,2),KMatrixSpace(K,4,4)!0)));
B:=VerticalJoin(B,HorizontalJoin(KMatrixSpace(K,4,8)!0,DnLat(4,2)));
x:=ChangeRing(1/2*Vector([K!1,1,1,1]),KK);
y:=ChangeRing(Vector([1,0,0,0]),KK);
B:=VerticalJoin(B,KMatrixSpace(K,1,12)![Eltseq(x+y) cat Eltseq(y) cat Eltseq(y)]);
B:=VerticalJoin(B,KMatrixSpace(K,1,12)![Eltseq(omega^2*x+y) cat Eltseq(omega*x+y) cat Eltseq(0*x)]);
B:=VerticalJoin(B,KMatrixSpace(K,1,12)![Eltseq(omega*x+y) cat Eltseq(0*x) cat Eltseq(omega^2*x+y)]);
lat1211bas:=ChangeRing(Matrix(TrueBasis(Module(Rows(ChangeRing(B,KK))))),K);
lat1211:=HermitianLattice(lat1211bas*HermitianTranspose(lat1211bas));

Lat12:=[lat121,lat122,lat123,lat124,lat125,lat126,lat127,lat128,lat129,lat1210,lat1211];

Unimod12:=Lat12 cat [DirectSum(II(1),x): x in Unimod11] cat [DirectSum(U6,U6)];