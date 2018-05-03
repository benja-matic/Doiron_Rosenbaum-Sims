#include("brent_riv.jl")
include("BRI.jl")
include("c://Users\\cohenbp\\Documents\\Neuroscience\\Balance_Circuits\\two_pools\\Analyze.jl")

using Distributions
srand(1234)

# ER1 = []
# ER2 = []
# IR1 = []
# IR2 = []
# TRE = []
# TRI = []
# wta = []
# for i in [2000, 4000, 8000, 16000, 32000]

N = 5000
half = div(N,2)
quar = div(half,2)
runtime = 10000
rt = runtime /1000.
h = 0.1
ntotal = Int64(runtime/h)

Jee = 12.5
Jie = 20.
Jei = -50.
Jii = -50.

vth = 20.
tau_m = 20.
tau_s = 2.
tau_a = 1000.
g_a = 5000.0/tau_a

p = 1/8.
k = 2500

BW = Brent_W(N, k, Jee, Jie, Jei, Jii)
CSR = sparse_rep(BW, N)

s1 = zeros(ntotal) .+ 2.4*h
s2 = zeros(ntotal) .+ 2.*h

# s1[30001:ntotal] = 2.7*h


@time t, r = Brent_Network_Euler_CSR_A(h, runtime, CSR, BW, N, s1, s2, vth, tau_m, tau_s, tau_a, g_a)

subplot(211)
title("Brent Network in Response to Various Inputs")
ylabel("Neuron #")
plot(t ./ 10000, r, ".", ms = 1.)
subplot(212)
plot(linspace(0,rt,ntotal), s1, alpha = .5)
plot(linspace(0,rt,ntotal), s2, alpha = .5)
plot(0, 0, "r.", alpha = 0.)
xlabel("Time")
ylabel("Input (mV/ms)")

# e_m = find(r .<= half);
# i_m = find(r .> half);
# ER = length(e_m)*(1/rt)*(1/half);
# IR = length(i_m)*(1/rt)*(1/half);
# te_pt = t[e_m];
# re_pt = r[e_m];
# ti_pt = t[i_m];
# ri_pt = r[i_m];
#
# e1 = find(te_pt .> quar)
# e2 = find(te_pt .<= quar)
# te1 = te_pt[e1]
# re1 = re_pt[e1]
# te2 = te_pt[e2]
# re2 = re_pt[e2]
#
# netd_bins = collect(1:netd_binsize:ntotal)
# nta = zeros(length(netd_bins)-1)
# ntb = zeros(length(netd_bins)-1)
#
# for j = 2:length(netd_bins)
#   tf = find(netd_bins[j-1] .<= te_pt .< netd_bins[j])
#   T = sum(re_pt[tf] .> quar)
#   B = sum(re_pt[tf] .<= quar)
#   nta[j-1] = T
#   ntb[j-1] = B
# end
#
# ntas  = 10000. .* (nta ./ quar) ./netd_binsize
# ntbs  = 10000. .* (ntb ./ quar) ./netd_binsize


# subplot(211)
# title("Strong Adaptation")
# plot(t, r, ".", ms = 1.)
# subplot(212)
# ylabel("Mean Firing Rate")
# plot(ntas)
# plot(ntbs)

#
# wta_ness, bias = score_analysis(re_pt, half)
# push!(wta, wta_ness)
#
# # eight = div(quar, 2)
# # plot(t/10000., r, "g.", ms = 1.)
# #
# # pools = ["E1", "E2", "I1", "I2"]
# # y = [eight, eight + quar, eight + half, eight + half + quar]
# # for i = 1:length(y)
# #              annotate(pools[i], xy = (-.4, y[i]))
# # end
# # ylabel("Neuron Index")
# # xlabel("Time (s)")
# #
# # title("Random Network w/ Asymmetric Inputs\ng_a=$(g_a)")
# # function s2ta(A, k, tau_s)
# #     return A*sqrt(k)*tau_s/1000.
# # end
#
# Wee=s2ta(abs(Jee), k, tau_s)
# Wei=s2ta(abs(Jei), k, tau_s)
# Wie=s2ta(abs(Jie), k, tau_s)
# Wii=s2ta(abs(Jii), k, tau_s)
#
# RE1_2x2, RI1_2x2 = theory_rates_2x2(Wee, Wie, Wei, Wii, max(s1, s2), max(s1, s2)-.1)
# E_R_top = [length(find(re_pt .== i))/rt for i=quar+1:half]
# E_R_bot = [length(find(re_pt .== i))/rt for i=1:quar]
# I_R_top = [length(find(ri_pt .== i))/rt for i=half+1:half+quar]
# I_R_bot = [length(find(ri_pt .== i))/rt for i=half+quar+1:N]
# #
# MER1 = mean(E_R_bot)#*(1/1000.)
# MER2 = mean(E_R_top)#*(1/1000.)
# MIR1 = mean(I_R_bot)#*(1/1000.)
# MIR2 = mean(I_R_top)#*(1/1000.)
#
# # push!(ER1, MER1)
# # push!(ER2, MER2)
# # push!(IR1, MIR1)
# # push!(IR2, MIR2)
# # push!(TRE, RE1_2x2)
# # push!(TRI, RI1_2x2)
#
# # end
# #
# plot([2000, 4000, 8000, 16000, 32000], ER1, ".", ms = 10., label = "Sim: ER1")
# plot([2000, 4000, 8000, 16000, 32000], ER2, ".", ms = 10., label = "Sim: ER2")
# plot([2000, 4000, 8000, 16000, 32000], IR1, ".", ms = 10., label = "Sim: IR1")
# plot([2000, 4000, 8000, 16000, 32000], IR2, ".", ms = 10., label = "Sim: IR2")
# plot([2000, 4000, 8000, 16000, 32000], TRE, ".", ms = 10., label = "2x2 theory: ER1")
# plot([2000, 4000, 8000, 16000, 32000], TRI, ".", ms = 10., label = "2x2 theory: IR1")
# xlabel("N")
# ylabel("Rate: Hz")
