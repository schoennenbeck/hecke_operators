IntertwiningOperator:=function(GenL, P)
    AutGenL := [AutomorphismGroup(L): L in GenL]
    SS:=Neighbours(GenL[1],P:AutoOrbits:=AutGenL[1])[1] meet GenL[1];
    PNB:=[[]: i in [1..#GenL]];
    SubGen:=[];
    for i in [1..#GenL] do
        M:=GenL[i];
        SM:=[S: S in MaximalSublattices(M,P: AutoOrbits:=AutGenL[i]) | IsLocallyIsometric(S,SS,P)];
        for S in SM do
            found:=false;
            for j in [1..#SubGen] do
                if IsIsometric(SubGen[j],S) then
                    PNB[i][j] +:= #AutGenL[i]/#(AutGenL[i] meet AutomorphismGroup(S));
                    found := true;
                    continue S;
                end if;
            end for;
            if not found then
                Append(~SubGen,S);
                Append(~PNB[i],#AutGenL[i]/#(AutGenL[i] meet AutomorphismGroup(S)));
                for k in [1..i-1] cat [i+1..#GenL] do
                    Append(~PNB[k],0);
                end for;
            end if;  
        end for;
    end for;
    return Matrix(PNB), SubGen
end function;

HeckeOperator:=function(L, P)
    GenL := GenusRepresentatives(L);
    AutGenL := [AutomorphismGroup(M): M in GenL]
    T21, SubGen := IntertwiningOperator(GenL, P)
    AutSubGen:=[AutomorphismGroup(x): x in SubGen;
    D0:=DiagonalMatrix([1/#x: x in AutGenL]);
    D1:=DiagonalMatrix([1/#x: x in AutSubGen]);
    T21:=ChangeRing(Matrix(T21),Rationals());
    T12:= D1^-1 * Transpose(T21) * D0;
    H := (MatrixRing(Rationals(), Nrows(T21))!T21*T12) - &+Eltseq(T21[1])
    return H, T21, T12, GenL, SubGen
end function;