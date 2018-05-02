#
include("brent_riv.jl")
include("c://Users\\cohenbp\\Documents\\Neuroscience\\Balance_Circuits\\two_pools\\Analyze.jl")

using Distributions
srand(1234)

N = 4000
half = N/2.
quar = half/2.
runtime = 30000
rt = runtime /1000.
h = 0.1

Jee = 12.5
Jie = 20.
Jei = -50.
Jii = -50.

vth = 20.
tau_m = 20.
tau_s = 2.
tau_a = 5000.
g_a = 3000.0/tau_a

BW = Brent_W(N, 1/8., Jee, Jie, Jei, Jii)
CSR = sparse_rep(BW, N)

s1 = 2.7
s2 = 2.4

t, r = Brent_Network_Euler_CSR_A(h, runtime, CSR, BW, N, s1, s2, vth, tau_m, tau_s, tau_a, g_a)

e_m = find(r .<= half);
i_m = find(r .> half);
ER = length(e_m)*(1/rt)*(1/half);
IR = length(i_m)*(1/rt)*(1/half);
te_pt = t[e_m];
re_pt = r[e_m];
ti_pt = t[i_m];
ri_pt = r[i_m];

wta_ness, bias = score_analysis(re_pt, half)
push!(wta, wta_ness)

eight = div(quar, 2)
plot(t/10000., r, "g.", ms = 1.)

pools = ["E1", "E2", "I1", "I2"]
y = [eight, eight + quar, eight + half, eight + half + quar]
for i = 1:length(y)
             annotate(pools[i], xy = (-.4, y[i]))
end
ylabel("Neuron Index")
xlabel("Time (s)")

title("Random Network w/ Asymmetric Inputs\ng_a=$(g_a)")
