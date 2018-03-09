#
include("brent_local.jl")
using Distributions
srand(1234)

N = 4000
half = N/2.
quar = half/2.
runtime = 10000
rt = runtime /1000.
h = 0.1
rot = round(Int64, runtime/h)
rot2 = rot*2

stdev = .5
tau_n = 200
mu = .2

R1 = rand(Normal(0, stdev), rot2)
R2 = rand(Normal(0, stdev), rot2)

s1 = OU_Model(R1, tau_n, 0.1)
s2 = OU_Model(R2, tau_n, 0.1)

s1 .+= mu
s2 .+= mu

s1 = s1[rot:end]
s2 = s2[rot:end]

BW = Brent_W(N, 1/8.)
CSR = sparse_rep(BW, N)

@time t1, r1 = Brent_Network_Local_Euler_CSR(h, runtime, CSR, BW, N, s1, 20., 20., 2.)
@time t2, r2 = Brent_Network_Local_Euler_CSR(h, runtime, CSR, BW, N, s2, 20., 20., 2.)

e_m1 = find(r1 .<= half)
i_ml = find(r1 .> half)

te1 = t1[e_m1]
re1 = r1[e_m1]

e_m2 = find(r2 .<= half)
i_m2 = find(r2 .> half)

te2 = t2[e_m2]
re2 = r2[e_m2]

r1T = length(find(re1 .> quar))
r1T = length(find(re2 .> quar))
