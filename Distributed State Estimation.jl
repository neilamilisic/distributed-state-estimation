
@everywhere using Distributed
@everywhere using SparseArrays
@everywhere using Printf
@everywhere using SuiteSparse
@everywhere using Random
@everywhere using StatsBase
@everywhere using LinearAlgebra
@everywhere using StatsBase


if nprocs()==1
addprocs(4)
end

@everywhere function matricaPovezanostiSistema(fromVec, toVec)
    n=max(maximum(fromVec), maximum(toVec))
    Y=spzeros(n,n)
    for i in eachindex(fromVec)
        Y[fromVec[i], toVec[i]] = 1
        Y[toVec[i], fromVec[i]] = 1
    end
return Y
end

@everywhere function pripada(n, S)
    for i in eachindex(S)
        if S[i]==n
        return true
        end
    end
    return false
end

@everywhere function matricaPovezanostiPodsistema(P, S)
    n=length(S)
    Y=spzeros(n,n)

    for j in 1:(length(P.colptr)-1)
        for i in P.colptr[j]:(P.colptr[j+1]-1)
            if P.nzval[i]==1 && pripada(P.rowval[i],S)==true && pripada(j,S) ==true
                a=findfirst(isequal(P.rowval[i]), S)
                b=findfirst(isequal(j), S)
                Y[a,b]=1

            end
        end

    end
return Y
end

@everywhere function interniJakobijan(mjerSnage, pseudoMjer, P)
        n=length(mjerSnage.nzval)+length(pseudoMjer.nzval)
        J=spzeros(n,size(P,1))
        k=length(pseudoMjer.nzval)
        for i in 1:length(pseudoMjer.nzval)
            J[i,findnz(pseudoMjer)[1][i]]=pseudoMjer.nzval[i]
        end
        br=k+1

        for j in 1:(length(mjerSnage.colptr)-1)
            for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                if j!=mjerSnage.rowval[i]
                    J[br,mjerSnage.rowval[i]]=1
                    J[br,j]=-1
                    br=br+1
                end
            end

        end

        for j in 1:(length(mjerSnage.colptr)-1)
            for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                if j==mjerSnage.rowval[i]
                    suma=0
                    for l in 1:(length(P.colptr)-1)
                        for r in P.colptr[l]:(P.colptr[l+1]-1)
                            if P.rowval[r]==j && P[P.rowval[r],l]!=0
                            J[br, l]=-1
                            suma=suma+1
                            end
                        end
                    end
                J[br,j]=suma
                br=br+1
                end

            end
        end
    return J
end

@everywhere function granicniJakobijan(mjerSnage, pseudoMjer, P, S, g)
        n=length(mjerSnage.nzval)+length(pseudoMjer.nzval)
        J=spzeros(n,length(S[g]))
        c=1

    for z in 1:length(S)

            for j in 1:(length(mjerSnage.colptr)-1)
                for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                    if j!=mjerSnage.rowval[i]
                        if pripada(mjerSnage.rowval[i], S[z])==true
                            if pripada(mjerSnage.rowval[i], S[g])==true
                                b=findfirst(isequal(mjerSnage.rowval[i]), S[g])
                                J[c,b]=1
                            end
                            if pripada(j, S[g])==true
                                b=findfirst(isequal(j), S[g])
                                J[c,b]=-1
                            end
                            c=c+1
                        end
                    end
                end
            end

            for j in 1:(length(mjerSnage.colptr)-1)
                for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                    if j==mjerSnage.rowval[i]
                        if pripada(mjerSnage.rowval[i], S[z])==true
                            if pripada(mjerSnage.rowval[i], S[g])==true
                                f=P[j,:]
                                w=length(f.nzval)
                                nz=findall(!iszero, f)

                                b=findfirst(isequal(mjerSnage.rowval[i]), S[g])
                                J[c,b]=w

                                for t in eachindex(nz)

                                    if pripada(nz[t], S[g])==true
                                        y=findfirst(isequal(nz[t]), S[g])
                                        J[c,y]=-1
                                    end
                                end
                            end

                            if pripada(mjerSnage.rowval[i], S[g])==false
                                f=P[j,:]
                                w=length(f.nzval)
                                nz=findall(!iszero, f)
                                for t in eachindex(nz)
                                    if pripada(nz[t], S[g])==true
                                        y=findfirst(isequal(nz[t]), S[g])
                                        J[c,y]=-1
                                    end
                                end
                            end
                        c=c+1
                        end
                    end
                end
            end

        ind=findall(!iszero, pseudoMjer)
        for i in 1:length(ind)
            if pripada(ind[i],S[z])==true
                if pripada(ind[i],S[g])==true
                    o=findfirst(isequal(ind[i]), S[g])
                    J[c,o]=1
                end
            c=c+1

            end
        end

    end
    return J
end

@everywhere function granicneKovarijanse(mjerSnage, pseudoMjer, S, g)
    n=length(mjerSnage.nzval)+length(pseudoMjer.nzval)
    R=spzeros(n,n)
    c=1
    sigma=0.01
for z in 1:length(S)

        for j in 1:(length(mjerSnage.colptr)-1)
            for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                if j!=mjerSnage.rowval[i]
                    if pripada(mjerSnage.rowval[i], S[z])==true
                        if pripada(mjerSnage.rowval[i], S[g])==true
                        R[c,c]=sigma*sigma
                        end
                        c=c+1
                    end
                end
            end
        end

        for j in 1:(length(mjerSnage.colptr)-1)
            for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                if j==mjerSnage.rowval[i]
                    if pripada(mjerSnage.rowval[i], S[z])==true
                        if pripada(mjerSnage.rowval[i], S[g])==true
                        R[c,c]=sigma*sigma
                        end
                    c=c+1
                    end
                end
            end
        end

    ind=findall(!iszero, pseudoMjer)
    for i in 1:length(ind)
        if pripada(ind[i],S[z])==true
            if pripada(ind[i], S[g])==true
            R[c,c]=-(sigma*sigma)
            end
        c=c+1
        end
    end

end
    return R
end

@everywhere function obzervabilnostSistema(H)
if rank(H)>=size(H,2)-1
        return true
        else return false
    end
end

@everywhere function saberiKolone(M)
    s=spzeros(size(M,1))
    for i in 1:size(M,2)
        s=M[:,i]+s
    end

    return s
end
@everywhere function brisanjePseudoMjerenja(M, mjerSnage, pseudoMjer, S)
    #brišu se psudo mjerenja
    c=1
        for z in 1:length(S)

                for j in 1:(length(mjerSnage.colptr)-1)
                    for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                        if j!=mjerSnage.rowval[i]
                            if pripada(mjerSnage.rowval[i], S[z])==true
                                c=c+1
                            end
                        end
                    end
                end

                for j in 1:(length(mjerSnage.colptr)-1)
                    for i in mjerSnage.colptr[j]:(mjerSnage.colptr[j+1]-1)
                        if j==mjerSnage.rowval[i]
                            if pripada(mjerSnage.rowval[i], S[z])==true
                            c=c+1
                            end
                        end
                    end
                end

            ind=findall(!iszero, pseudoMjer)
            for i in 1:length(ind)
                if pripada(ind[i],S[z])==true
                    for j in 1:size(M,2)
                    M[c,j]=0
                    end
                c=c+1
                end
            end
        end
        dropzeros!(M)
        return M
end

@everywhere function r_i(mjerSnage, pseudoMjer, delta, z, J)
    h=zeros(length(z))
        r=z-J*delta*(pi/180)
        return r
end

@everywhere function r_g(mjerSnage, pseudoMjer, delta, z, J1_g, J2_g, J3_g, J4_g)
    J=hcat(J1_g, J2_g, J3_g, J4_g)
        r=z-J*delta*(pi/180)
        return r
end

@everywhere function diagsqrt(M)
    v=spzeros(size(M,1),size(M,2))
    for i in 1:size(M,1)
        v[i][i]=sqrt(M[i][i])
    end
    return v
end

@everywhere function max_element(r_n)
    max_e=r_n[1][1]
    n=1
    m=1
    for i in 1:length(r_n)
        for j in 1:length(r_n[i])
            if max_e<r_n[i][j]
            maxi_e=r_n[i][j]
            n=i
            m=j
            end
        end
    end
    return n, m
end
#-----------------------Ulazni podaci za IEEE14-------------------------------

z1=[0.0000, 0.0869, 0.3532, 0.0663, 0.2402]
z2=[0.0000, -0.0243, 0.0356, -0.0002]
z3=[0.0000, 0.0099, 0.0148, 0.0164, 0.0016, -0.0133]
z4=[0.0000, 0.0028, 0.0192]
zc1=[-0.0801]
zc2=[-0.0445, 0.0630, 0.0274, -0.1594, 0.0000]
zc3=[0.0154, -0.2026, 0.0000]
zc4=[-0.0054, -0.0346, 0.0000]
zc=[-0.0801, -0.0445, 0.0630, 0.0274, -0.1594, 0.0000, 0.0154, -0.2026, 0.0000, -0.0054, -0.0346, 0.0000]
#Podsistemi
S1=[1,2,5]
S2=[3,4,7,8]
S3=[6,11,12,13]
S4=[9,10,14]

S=[S1,S2,S3,S4]

#Matrica povezanosti sistema
fromVec=[1,1,2,2,2,3,4,4,4,5,6,6,6,7,7,9,9,10,12,13]
toVec=[2,5,3,4,5,4,5,7,9,6,11,12,13,8,9,10,14,11,13,14]
P=@spawnat 1 matricaPovezanostiSistema(fetch(fromVec), fetch(toVec))
#Podsistem1
P1=@spawnat 2 matricaPovezanostiPodsistema(fetch(P),fetch(S1))
#Podsistem2
P2=@spawnat 3 matricaPovezanostiPodsistema(fetch(P), fetch(S2))
#Podsistem3
P3=@spawnat 4 matricaPovezanostiPodsistema(fetch(P), fetch(S3))
#Podsistem4
P4=@spawnat 5 matricaPovezanostiPodsistema(fetch(P), fetch(S4))


#= Matrice mjerenja snaga INTERNA
* 1 na dijagonali - mjerenje injektirane aktivne snage
* 1 van dijagonale - mjerenje toka aktivne snage =#
mjerSnage1= @spawnat 2 sparse([1, 1, 1, 2], [1, 2, 3, 3], [1, 1, 1, 1])
mjerSnage2= @spawnat 3 sparse([1, 2, 3], [2, 3, 4], [1, 1, 1])
mjerSnage3= @spawnat 4 sparse([1, 1, 1, 3, 3], [2, 3, 4, 3, 4], [1, 1, 1, 1, 1])
mjerSnage4= @spawnat 5 sparse([1, 1], [2, 3], [1, 1])

#= Vektori pseudo mjerenja INTERNI
* 1 - pseudo mjerenje na odgovarajućem čvoru =#
pseudoMjer1=@spawnat 2 sparsevec([1],[1])
pseudoMjer2=@spawnat 3 sparsevec([1],[1])
pseudoMjer3=@spawnat 4 sparsevec([1],[1])
pseudoMjer4=@spawnat 5 sparsevec([1],[1])

#= Matrice mjerenja snaga GRANIČNA
* 1 na dijagonali - mjerenje injektirane aktivne snage
* 1 van dijagonale - mjerenje toka aktivne snage =#
mjerSnage_g= @spawnat 1 sparse([5,4,4,7,3,13,13,10,14], [5,5,9,9,3,14,13,11,14], [1,1,1,1,1,1,1,1,1])
#= Vektori pseudo mjerenja GRANIČNI
* 1 - pseudo mjerenje na odgovarajućem čvoru =#
pseudoMjer_g=@spawnat 1 sparsevec([3,6,9],[1,1,1])
#-----------------------------Jakobijani internih mjerenja-------------------------------------------------------------
J1=@spawnat 2 interniJakobijan(fetch(mjerSnage1), fetch(pseudoMjer1), fetch(P1))
J2=@spawnat 3 interniJakobijan(fetch(mjerSnage2), fetch(pseudoMjer2), fetch(P2))
J3=@spawnat 4 interniJakobijan(fetch(mjerSnage3), fetch(pseudoMjer3), fetch(P3))
J4=@spawnat 5 interniJakobijan(fetch(mjerSnage4), fetch(pseudoMjer4), fetch(P4))
#----------------------------Jakobijani graničnih mjerenja-------------------------------------------------------------
J1_g=@spawnat 2 granicniJakobijan(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(P), fetch(S),fetch(1))
J2_g=@spawnat 3 granicniJakobijan(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(P), fetch(S),fetch(2))
J3_g=@spawnat 4 granicniJakobijan(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(P), fetch(S),fetch(3))
J4_g=@spawnat 5 granicniJakobijan(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(P), fetch(S),fetch(4))
#----------------------------Granične kovarijansne matrice-------------------------------------------------------------
R1_g=@spawnat 2 granicneKovarijanse(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(S), 1)
R2_g=@spawnat 3 granicneKovarijanse(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(S), 2)
R3_g=@spawnat 4 granicneKovarijanse(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(S), 3)
R4_g=@spawnat 5 granicneKovarijanse(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(S), 4)



#---------------------------Lokalne kovarijansne matrice---------------------------------------------------------------
sigma=0.01
rl=sigma*sigma

R1=@spawnat 2 sparse(rl*I, maximum(fetch(findnz(fetch(J1))[1])), maximum(fetch(findnz(fetch(J1))[1])))
R2=@spawnat 3 sparse(rl*I, maximum(fetch(findnz(fetch(J2))[1])), maximum(fetch(findnz(fetch(J2))[1])))
R3=@spawnat 4 sparse(rl*I, maximum(fetch(findnz(fetch(J3))[1])), maximum(fetch(findnz(fetch(J3))[1])))
R4=@spawnat 5 sparse(rl*I, maximum(fetch(findnz(fetch(J4))[1])), maximum(fetch(findnz(fetch(J4))[1])))
#---------------------------------Gain matrice-------------------------------------------------------------------------
G1=@spawnat 2 transpose(fetch(J1))*(fetch(R1)\one(fetch(R1)))*fetch(J1)
G2=@spawnat 3 transpose(fetch(J2))*(fetch(R2)\one(fetch(R2)))*fetch(J2)
G3=@spawnat 4 transpose(fetch(J3))*(fetch(R3)\one(fetch(R3)))*fetch(J3)
G4=@spawnat 5 transpose(fetch(J4))*(fetch(R4)\one(fetch(R4)))*fetch(J4)
#-------------------Provjera obzervalbilnosti podsistema----------------------------------------------------------------
if obzervabilnostSistema(fetch(G1))==false
    println("Podsistem 1 nije obzervabilan!")
    return
end
if obzervabilnostSistema(fetch(G2))==false
    println("Podsistem 2 nije obzervabilan!")
    return
end
if obzervabilnostSistema(fetch(G3))==false
    println("Podsistem 3 nije obzervabilan!")
    return
end
if obzervabilnostSistema(fetch(G4))==false
    println("Podsistem 4 nije obzervabilan!")
    return
end
#------------------------------Lokalna estimacija (fazni uglovi pretvoreni iz radijana u stepene)---------------------------
theta1=@spawnat 2 (fetch(G1)\Matrix(one(fetch(G1))))*transpose(fetch(J1))*(fetch(R1)\one(fetch(R1)))*fetch(z1)*(180/pi)
theta2=@spawnat 3 (fetch(G2)\Matrix(one(fetch(G2))))*transpose(fetch(J2))*(fetch(R2)\one(fetch(R2)))*fetch(z2)*(180/pi)
theta3=@spawnat 4 (fetch(G3)\Matrix(one(fetch(G3))))*transpose(fetch(J3))*(fetch(R3)\one(fetch(R3)))*fetch(z3)*(180/pi)
theta4=@spawnat 5 (fetch(G4)\Matrix(one(fetch(G4))))*transpose(fetch(J4))*(fetch(R4)\one(fetch(R4)))*fetch(z4)*(180/pi)
#-------------------Provjera obzervabilnosti sistema-------------------------------------------------------------------
W1=@spawnat 2 saberiKolone(fetch(J1_g))
W2=@spawnat 3 saberiKolone(fetch(J2_g))
W3=@spawnat 4 saberiKolone(fetch(J3_g))
W4=@spawnat 5 saberiKolone(fetch(J4_g))
W=hcat(fetch(W1),fetch(W2),fetch(W3),fetch(W4))
Wc=brisanjePseudoMjerenja(fetch(W), fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(S))
if obzervabilnostSistema(fetch(W))==false
    println("Sistem nije obzervabilan!")
    return
end
#---------------------------G_ci matrice-------------------------------------------------------------------------------
G1_g=@spawnat 2 fetch(R1_g)+fetch(J1_g)*(fetch(G1)\Matrix(one(fetch(G1))))*transpose(fetch(J1_g))
G2_g=@spawnat 3 fetch(R2_g)+fetch(J2_g)*(fetch(G2)\Matrix(one(fetch(G2))))*transpose(fetch(J2_g))
G3_g=@spawnat 4 fetch(R3_g)+fetch(J3_g)*(fetch(G3)\Matrix(one(fetch(G3))))*transpose(fetch(J3_g))
G4_g=@spawnat 5 fetch(R4_g)+fetch(J4_g)*(fetch(G4)\Matrix(one(fetch(G4))))*transpose(fetch(J4_g))
#---------------------Podaci se proslijeđuju do koordinatora gdje se formira Gc-----------------------------------------
G_g=@spawnat 1 fetch(G1_g)+fetch(G2_g)+fetch(G3_g)+fetch(G4_g)
#------------------Vektor Lagrange-ovih multiplikatora------------------------------------------------------------------
lambda= @spawnat 1 (fetch(G_g)\Matrix(one(fetch(G_g))))*(fetch(zc)*(180/pi)-(fetch(J1_g)*fetch(theta1)+fetch(J2_g)*fetch(theta2)+fetch(J3_g)*fetch(theta3)+fetch(J4_g)*fetch(theta4)))
#------------------Lokalna estimacija na osnovu graničnih mjerenja u stepenima------------------------------------------------------------------------------------------------
fi1=@spawnat 2 (fetch(G1)\Matrix(one(fetch(G1))))*transpose(fetch(J1_g))*fetch(lambda)
fi2=@spawnat 3 (fetch(G2)\Matrix(one(fetch(G2))))*transpose(fetch(J2_g))*fetch(lambda)
fi3=@spawnat 4 (fetch(G3)\Matrix(one(fetch(G3))))*transpose(fetch(J3_g))*fetch(lambda)
fi4=@spawnat 5 (fetch(G4)\Matrix(one(fetch(G4))))*transpose(fetch(J4_g))*fetch(lambda)
#------------------Estimacija stanja cjelokupnog sistema------------------------------------------------------------------------------------------------------------------------
delta1=fetch(@spawnat 1 (fetch(theta1)+fetch(fi1)))
delta2=fetch(@spawnat 1 (fetch(theta2)+fetch(fi2)))
delta3=fetch(@spawnat 1 (fetch(theta3)+fetch(fi3)))
delta4=fetch(@spawnat 1 (fetch(theta4)+fetch(fi4)))
#-------------------Bad data detection------------------------------------------------------------------------------------------------------------------------------------------
r_1=@spawnat 2 r_i(fetch(mjerSnage1), fetch(pseudoMjer1), fetch(delta1), fetch(z1), fetch(J1))
r_2=@spawnat 3 r_i(fetch(mjerSnage2), fetch(pseudoMjer2), fetch(delta2), fetch(z2), fetch(J2))
r_3=@spawnat 4 r_i(fetch(mjerSnage3), fetch(pseudoMjer3), fetch(delta3), fetch(z3), fetch(J3))
r_4=@spawnat 5 r_i(fetch(mjerSnage4), fetch(pseudoMjer4), fetch(delta4), fetch(z4), fetch(J4))

delta_g=vcat(fetch(delta1),fetch(delta2),fetch(delta3),fetch(delta4))

r_c=@spawnat 1 r_g(fetch(mjerSnage_g), fetch(pseudoMjer_g), fetch(delta_g), fetch(zc),fetch(J1_g), fetch(J2_g), fetch(J3_g), fetch(J4_g))

P_1a=@spawnat 2 fetch(R1)-fetch(J1)*(fetch(G1)\Matrix(one(fetch(G1))))*transpose(fetch(J1))
P_1b=@spawnat 2 fetch(J1)*(fetch(G1)\Matrix(one(fetch(G1))))*(transpose(fetch(J1_g))*(fetch(G_g)\Matrix(one(fetch(G_g))))*fetch(J1_g))*(fetch(G1)\Matrix(one(fetch(G1))))*transpose(fetch(J1))
P_1=@spawnat 2 fetch(P_1a)+fetch(P_1b)

P_2a=@spawnat 3 fetch(R2)-fetch(J2)*(fetch(G2)\Matrix(one(fetch(G2))))*transpose(fetch(J2))
P_2b=@spawnat 3 fetch(J2)*(fetch(G2)\Matrix(one(fetch(G2))))*(transpose(fetch(J2_g))*(fetch(G_g)\Matrix(one(fetch(G_g))))*fetch(J2_g))*(fetch(G2)\Matrix(one(fetch(G2))))*transpose(fetch(J2))
P_2=@spawnat 3 fetch(P_2a)+fetch(P_2b)

P_3a=@spawnat 4 fetch(R3)-fetch(J3)*(fetch(G3)\Matrix(one(fetch(G3))))*transpose(fetch(J3))
P_3b=@spawnat 4 fetch(J3)*(fetch(G3)\Matrix(one(fetch(G3))))*(transpose(fetch(J3_g))*(fetch(G_g)\Matrix(one(fetch(G_g))))*fetch(J3_g))*(fetch(G3)\Matrix(one(fetch(G3))))*transpose(fetch(J3))
P_3=@spawnat 4 fetch(P_3a)+fetch(P_3b)

P_4a=@spawnat 5 fetch(R4)-fetch(J4)*(fetch(G4)\Matrix(one(fetch(G4))))*transpose(fetch(J4))
P_4b=@spawnat 5 fetch(J4)*(fetch(G4)\Matrix(one(fetch(G4))))*(transpose(fetch(J4_g))*(fetch(G_g)\Matrix(one(fetch(G_g))))*fetch(J4_g))*(fetch(G4)\Matrix(one(fetch(G4))))*transpose(fetch(J4))
P_4=@spawnat 5 fetch(P_4a)+fetch(P_4b)

R_c=@spawnat 1 fetch(R1_g)+fetch(R2_g)+fetch(R3_g)+fetch(R4_g)

P_c=@spawnat 1 fetch(R_c)*(fetch(G_g)\Matrix(one(fetch(G_g))))*fetch(R_c)

r1_n=fetch(@spawnat 2 (sqrt(complex(Diagonal(fetch(P_1))))\Matrix(one(sqrt(complex(Diagonal(fetch(P_1)))))))*fetch(r_1))
r2_n=fetch(@spawnat 3 (sqrt(complex(Diagonal(fetch(P_2))))\Matrix(one(sqrt(complex(Diagonal(fetch(P_2)))))))*fetch(r_2))
r3_n=fetch(@spawnat 4 (sqrt(complex(Diagonal(fetch(P_3))))\Matrix(one(sqrt(complex(Diagonal(fetch(P_3)))))))*fetch(r_3))
r4_n=fetch(@spawnat 5 (sqrt(complex(Diagonal(fetch(P_4))))\Matrix(one(sqrt(complex(Diagonal(fetch(P_4)))))))*fetch(r_4))

rc_n=fetch(@spawnat 1 (sqrt(complex(Diagonal(fetch(P_c))))\Matrix(one(sqrt(complex(Diagonal(fetch(P_c)))))))*fetch(r_c))

r_n=[real(fetch(r1_n)), real(fetch(r2_n)), real(fetch(r3_n)), real(fetch(r4_n))]

n,m=max_element(fetch(r_n))
prrrr=findnz(fetch(mjerSnage1))

III=findnz(fetch(mjerSnage_g))
#=

if fetch(maxi)>3
    if n==1
        if m<=length(fetch(pseudoMjer1))
            fetch(pseudoMjer1).nzval[m]=0
            dropzeros!(fetch(pseudoMjer1))
        end
        if m>length(fetch(pseudoMjer1))
            fetch(mjerSnage1).nzval[m-length(fetch(pseudoMjer1)]=0
            dropzeros!(fetch(mjerSnage1))
        end
    end

    if n==2
        if m<=length(fetch(pseudoMjer2))
            fetch(pseudoMjer2).nzval[m]=0
            dropzeros!(fetch(pseudoMjer2))
        end
        if m>length(fetch(pseudoMjer2))
            fetch(mjerSnage2).nzval[m-length(fetch(pseudoMjer2)]=0
            dropzeros!(fetch(mjerSnage2))
        end
    end

    if n==3
        if m<=length(fetch(pseudoMjer3))
            fetch(pseudoMjer3).nzval[m]=0
            dropzeros!(fetch(pseudoMjer3))
        end
        if m>length(fetch(pseudoMjer3))
            fetch(mjerSnage3).nzval[m-length(fetch(pseudoMjer3)]=0
            dropzeros!(fetch(mjerSnage3))
        end
    end

    if n==4
        if m<=length(fetch(pseudoMjer4))
            fetch(pseudoMjer4).nzval[m]=0
            dropzeros!(fetch(pseudoMjer4))
        end
        if m>length(fetch(pseudoMjer4))
            fetch(mjerSnage4).nzval[m-length(fetch(pseudoMjer4)]=0
            dropzeros!(fetch(mjerSnage4))
        end
    end



end

plplp=findmax(fetch(r_n))

#max1=fetch(findmax(real(r1_n)))
#max2=fetch(findmax(real(r2_n)))
#max3=fetch(findmax(real(r3_n)))
#max4=fetch(findmax(real(r4_n)))
#maxc=fetch(findmax(real(rc_n)))    =#
