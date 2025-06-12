#####Calculate the Lyapunov exponents of the LaNMM
#####Check https://juliadynamics.github.io/ChaosTools.jl/dev/lyapunovs/
#####Adapted from Clusella
using ChaosTools, StaticArrays, UnPack, NPZ, PyPlot, DynamicalSystems, OrdinaryDiffEq
r0 = 2.5
va = 6.0
vb = 1.0
ρ = 0.56
function Sigm(v, vr) 
    2.0 * r0 / (1 + exp.(ρ * (vr - v)))
end 

function dSigm(v,vr)
    dsigm =2*r0*ρ*exp.(ρ*(vr-v))/((1+exp.(ρ*(vr-v)))^2)
end

A_AMPA = 3.25
A_GABAs = -22.0
A_GABAf = -30.0
a_AMPA = 100.0
a_GABAs = 50.0
a_GABAf = 220.0

C = [108.0, 33.75, 0.0, 135.0, 33.75, 70.0, 550.0, 0.0, 200.0, 100.0, 80.0, 200.0, 30.0]
A = [A_AMPA, A_AMPA, A_GABAs, A_AMPA, A_GABAf]
b = [a_AMPA, a_AMPA, a_GABAs, a_AMPA, a_GABAf]
w = [va, va, va, vb, va]

function LaNMM(u, p, t)
    p1= p[1];p2=p[2]
    y0,y1,y2,y3,y4,y5,y6,y7,y8,y9 = u

    f = zeros(Float64, 10)

    f[1]= y5
    f[2]= y6
    f[3]= y7
    f[4]= y8
    f[5]= y9
    
    f[6]= A[1] * b[1] * Sigm(C[1] * y1 + C[2] * y2 + C[11] * y3 + A[1] * p1 / b[1], w[1]) - 2 * b[1] * y5 - b[1]^2 * y0
    f[7]= A[2] * b[2] * Sigm(C[4] * y0, w[2]) - 2 * b[2] * y6 - b[2]^2 * y1
    f[8]= A[3] * b[3] * Sigm(C[5] * y0, w[3]) - 2 * b[3] * y7 - b[3]^2 * y2
    f[9]= A[4] * b[4] * Sigm(C[6] * y3 +  C[7] * y4 + C[12] * y0 + A[4] * p2 / b[4] , w[4]) - 2 * b[4] * y8 - b[4]^2 * y3
    f[10]=A[5] * b[5] * Sigm(C[9] * y3 + C[10] * y4 + C[13] * y0,w[5]) - 2 * b[5] * y9 - b[5]^2 * y4

    return SVector{10}(f)
end


function tandy_LaNMM(u, p, t)
    p1 = p[1]; p2 = p[2]
    y0,y1,y2,y3,y4,y5,y6,y7,y8,y9 = u
    jac = zeros(Float64, 10, 10)  # Initialize an array of SVector rows
    jac[1,6]=1;
    jac[2,7]=1;
    jac[3,8]=1;
    jac[4,9]=1;
    jac[5,10]=1;

    jac[6,6]=-2*b[1];
    jac[7,7]=-2*b[2];
    jac[8,8]=-2*b[3];
    jac[9,9] = -2*b[4];
    jac[10,10] = -2*b[5];

    jac[6,1]=-b[1]^2;
    jac[6,2]= A[1]*b[1]*C[1] *dSigm(C[1] * y1 + C[2] * y2 + C[11] * y3 + A[1] * p1 / b[1],w[1]);
    jac[6,3]= A[1]*b[1]*C[2] *dSigm(C[1] * y1 + C[2] * y2 + C[11] * y3 + A[1] * p1 / b[1],w[1]);
    jac[6,4]= A[1]*b[1]*C[11]*dSigm(C[1] * y1 + C[2] * y2 + C[11] * y3 + A[1] * p1 / b[1],w[1]);

    jac[7,1]= A[2]*b[2]*C[4] *dSigm(C[4] * y0,w[2]);
    jac[7,2]=-b[2]^2;

    jac[8,1]= A[3]*b[3]*C[5] *dSigm(C[5] * y0,w[3]);
    jac[8,3]=-b[3]^2;

    jac[9,1] = A[4]*b[4]*dSigm(C[6] * y3 +  C[7] * y4 + C[12] * y0 + A[4] * p2 / b[4] ,w[4])*C[12]; 
    jac[9,4] = A[4]*b[4]*dSigm(C[6] * y3 +  C[7] * y4 + C[12] * y0 + A[4] * p2 / b[4] ,w[4])*C[6] - b[4]^2;
    jac[9,5] = A[4]*b[4]*dSigm(C[6] * y3 +  C[7] * y4 + C[12] * y0 + A[4] * p2 / b[4] ,w[4])*C[7]; 

    jac[10,1] = A[5]*b[5]*dSigm(C[9] * y3 + C[10] * y4 + C[13] * y0,w[5])*C[13];
    jac[10,4] = A[5]*b[5]*dSigm(C[9] * y3 + C[10] * y4 + C[13] * y0,w[5])*C[9];
    jac[10,5] = A[5]*b[5]*dSigm(C[9] * y3 + C[10] * y4 + C[13] * y0,w[5])*C[10] - b[5]^2;
    return SMatrix{10,10}(jac)
end


u0 = fill(0.05,10) 
tspan = (0.0, 100.0)

function LaNMM_rule(u, p, t)
    p1= p[1];p2=p[2]
    y0,y1,y2,y3,y4,y5,y6,y7,y8,y9 = u
    f = zeros(eltype(u), 10)
    f[1]= y5
    f[2]= y6
    f[3]= y7
    f[4]= y8
    f[5]= y9
    f[6]= A[1] * b[1] * Sigm(C[1] * y1 + C[2] * y2 + C[11] * y3 + A[1] * p1 / b[1], w[1]) - 2 * b[1] * y5 - b[1]^2 * y0
    f[7]= A[2] * b[2] * Sigm(C[4] * y0, w[2]) - 2 * b[2] * y6 - b[2]^2 * y1
    f[8]= A[3] * b[3] * Sigm(C[5] * y0, w[3]) - 2 * b[3] * y7 - b[3]^2 * y2
    f[9]= A[4] * b[4] * Sigm(C[6] * y3 +  C[7] * y4 + C[12] * y0 + A[4] * p2 / b[4] , w[4]) - 2 * b[4] * y8 - b[4]^2 * y3
    f[10]=A[5] * b[5] * Sigm(C[9] * y3 + C[10] * y4 + C[13] * y0,w[5]) - 2 * b[5] * y9 - b[5]^2 * y4
    return SVector{10}(f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10])
end


parg1 = parse(Float64,ARGS[1])
parg2 = parse(Float64,ARGS[2])

pp= [parg1,parg2]
nmm = CoupledODEs(LaNMM, u0, pp,diffeq = (alg =RK4(),dt=1e-3,adaptive=false))
#diffeq = (alg = RK4(), dtmax = 1e-3,dtmin=1e-5,dt=1e-3))

tandy = TangentDynamicalSystem(nmm; J = tandy_LaNMM)
#lanmm = CoupledODEs(LaNMM_rule, u0, pp,diffeq = (alg = RK4(), dtmax = 1e-3,dtmin=1e-5,dt=1e-4))

ly =lyapunovspectrum(tandy, 2000,Δt=0.5,Ttr=2000)[1:2]
println(ly)
