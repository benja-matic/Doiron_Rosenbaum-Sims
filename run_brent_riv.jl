#
include("brent_riv.jl")
include("c://Users\\cohenbp\\Documents\\Neuroscience\\Balance_Circuits\\two_pools\\Analyze.jl")

using Distributions
srand(1234)

N = 4000
half = div(N,2)
quar = div(half,2)
runtime = 10000
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
g_a = 0.

p = 1/8.
k = 500

BW = Brent_W(N, k, Jee, Jie, Jei, Jii)
CSR = sparse_rep(BW, N)

s1 = 2.7
s2 = 2.4

# t, r = Brent_Network_Euler_CSR_A(h, runtime, CSR, BW, N, s1, s2, vth, tau_m, tau_s, tau_a, g_a)
#
# e_m = find(r .<= half);
# i_m = find(r .> half);
# ER = length(e_m)*(1/rt)*(1/half);
# IR = length(i_m)*(1/rt)*(1/half);
# te_pt = t[e_m];
# re_pt = r[e_m];
# ti_pt = t[i_m];
# ri_pt = r[i_m];
#
# wta_ness, bias = score_analysis(re_pt, half)
# push!(wta, wta_ness)
#
# eight = div(quar, 2)
# plot(t/10000., r, "g.", ms = 1.)
#
# pools = ["E1", "E2", "I1", "I2"]
# y = [eight, eight + quar, eight + half, eight + half + quar]
# for i = 1:length(y)
#              annotate(pools[i], xy = (-.4, y[i]))
# end
# ylabel("Neuron Index")
# xlabel("Time (s)")
#
# title("Random Network w/ Asymmetric Inputs\ng_a=$(g_a)")
# function s2ta(A, k, tau_s)
#     return A*sqrt(k)*tau_s/1000.
# end
#
# Wee=s2ta(Jee, round(Int64, N*p), 2)
# Wei=s2ta(Jei, round(Int64, N*p), 2)
# Wie=s2ta(Jie, round(Int64, N*p), 2)
# Wii=s2ta(Jii, round(Int64, N*p), 2)
#
# RE1_2x2, RI1_2x2 = theory_rates_2x2(Wee, Wie, Wei, Wii, max(s1, s2), max(s1, s2)-.1)
# E_R_top = [length(find(re_pt .== i))/rt for i=quar+1:half]
# E_R_bot = [length(find(re_pt .== i))/rt for i=1:quar]
# I_R_top = [length(find(ri_pt .== i))/rt for i=half+1:half+quar]
# I_R_bot = [length(find(ri_pt .== i))/rt for i=half+quar+1:N]
#
# MER1 = mean(E_R_bot)#*(1/1000.)
# MER2 = mean(E_R_top)#*(1/1000.)
# MIR1 = mean(I_R_bot)#*(1/1000.)
# MIR2 = mean(I_R_top)#*(1/1000.)
