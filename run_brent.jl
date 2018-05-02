#run brent small
include("brent.jl")
srand(1234)
N = 5000
half = round(Int64, N/2)
quar = round(Int64, N/4)

runtime = 10000
rt = runtime /1000.
h = 0.1
rot = round(Int64, runtime/h)
rot2 = rot*2
netd_binsize = 150/h
cbinsize = 1000.

BW = Brent_W(N, 1/8.)
CSR = sparse_rep(BW, N)

scan = linspace(.2, .8, 20)
MD = []
alt = []

stdev = .1*5
tau_n = 100*10
mu = 0.2

# R1 = rand(Normal(0, stdev), rot2)
# R2 = rand(Normal(0, stdev), rot2)
#
# s1 = OU_Model(R1, tau_n, 0.1)
# s2 = OU_Model(R2, tau_n, 0.1)
#
# s1 .+= mu
# s2 .+= mu
#
# s1 = s1[rot:end]
# s2 = s2[rot:end]
#
# STD = mean([std(s1), std(s2)])
# stdev/(sqrt(2)*sqrt(tau_n))
#
# A = s1./(s1 .+ s2)

s1 = zeros(rot) .+ .2#.35
s2 = zeros(rot) .+ .3#.35


#stdev_resulting = stdev/(sqrt(2)*sqrt(tau_n))

@time t, r = Brent_Network_Euler_CSR(h, runtime, CSR, BW, N, s1, s2, 20., 20., 2.)

e_m = find(r .<= half)
i_m = find(r .> half)


ER = length(e_m)*(1/rt)*(1/2500.)
IR = length(i_m)*(1/rt)*(1/2500.)

#
te_pt = t[e_m]
re_pt = r[e_m]
ti_pt = t[i_m]
ri_pt = r[i_m]
#
E_neurons = Neuron_finder(re_pt, 5, 5)
I_neurons = Neuron_finder(ri_pt, 5, 5)
E_Rates = [length(find(r .== i))/rt for i=1:half]
I_Rates = [length(find(r .== i))/rt for i=half+1:N]
CVE = CV_ISI_ALLTIME(E_neurons, te_pt, re_pt)
CVI = CV_ISI_ALLTIME(I_neurons, ti_pt, ri_pt)
RSC_TOP = rand_pair_cor(cbinsize, te_pt, re_pt, E_neurons, 1000)

# ER = []
# IR = []
# ECV = []
# ERC = []
# ICV = []

# push!(ER, mean(E_Rates))
# push!(IR, mean(I_Rates))
# push!(ECV, mean(CVE))
# push!(ICV, mean(CVI))
# push!(ERC, mean(RSC_TOP))

# ntd, nts = nt_diff(te_pt, re_pt, rot, quar, 25/h)
# s = ntd ./ nts

# push!(sm, mean(s))
# push!(ss, std(s))
# include("c://Users\\cohenbp\\Documents\\Neuroscience\\Balance_Circuits\\two_pools\\Analyze.jl")

# ntd, nts = nt_diff_H(te_pt, re_pt, rot, quar, netd_binsize)
# s = ntd ./ nts #signal for dominances
# flags, times = WLD_01(s, -.333, .333)
# tx = times .* netd_binsize
#
# d = convert(Array{Float64}, diff(tx ./ (1000. / h)))
# cvd = cv(d) ###Raw estimate of CVD, likely to include very rapid switches which should really be smoothed out
#
# LP = .3
#
# dx = []
# for i in d
#   if i > LP
#       push!(dx, i)
#   end
# end
# dx = convert(Array{Float64}, dx)
# cvdlp = cv(dx) ###Low-pass filter measure of CVD
#
# ###Use this code for spiking statistics, which includes only dominance and suppression times; mixed percept time is absorbed into either of those
# t2, f2 = splice_reversions(flags, times) ###Another way to get rid of rapid switches that aren't really there
# fw = find(f2 .== "win")
# fl = find(f2 .== "lose")
# tx2 = t2 .* netd_binsize
# d2 = convert(Array{Float64}, diff(tx2))
# cvd2 = cv(d2) ### Another estimate of CVD

# tl = .3
# th = .7
#
# flags, times = WLD_01(s, tl, th)
# top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)
# m = (rot - tnmz)/(length(top) + length(bot))
# MDT = tdom/length(top)
# MDB = bdom/length(bot)
# ad = []
# for i in eachindex(top)
#     push!(ad, top[i][2]-top[i][1])
# end
# for i in eachindex(bot)
#     push!(ad, bot[i][2]-bot[i][1])
# end
# ad = convert(Array{Float64}, ad)/10000.
# adx = find(ad .> 0.2)
# adz = ad[adx]
# cvd = cv(ad)
# println("##RESULT $(MDT), $(MDB), $(m), $(cvd), $(length(times))")
# push!(alt, length(times))
# push!(MD, mean(ad))
#
# plt[:hist](adz, 20)


#
#
# if flags == ["lose", "end"]
#     println("##RESULT $(0), $(rot), $(rot), $(0)")
# elseif flags == ["win", "end"]
#     println("##RESULT $(rot), $(0), $(rot), $(0)")
# elseif flags == ["draw", "end"]
#     println("##RESULT $(0), $(0), $(0), $(0)")
# else
#     t2, f2 = splice_reversions(flags, times)
#     top, tdom, bot, bdom, nmz, tnmz = splice_flags(flags, times)
#
#     # d = convert(Array{Float64}, diff(250*t2))
#     # m = mean(d)
#     # cvd = cv(d)
#     # fw = find(f2.=="win")
#     # fl = find(f2.=="lose")
#     # MDT = mean(d[fw])
#     # MDB = mean(d[fl])
#     m = (rot - tnmz)/(length(top) + length(bot))
#     MDT = tdom/length(top)
#     MDB = bdom/length(bot)
#     ad = []
#     for i in eachindex(top)
#         push!(ad, top[i][2]-top[i][1])
#     end
#     for i in eachindex(bot)
#         push!(ad, bot[i][2]-bot[i][1])
#     end
#     ad = convert(Array{Float64}, ad)/10000.
#     cvd = cv(ad)
#
#     println("##RESULT $(MDT), $(MDB), $(m), $(cvd)")
# end

# subplot(311)
# ax1 = gca()
# title("Simulation of Brent's Network")
# plot(te_pt./10000, re_pt, "g.", ms = 1.)
# plot(t2*25/1000, fill(quar, length(t2)), "k.", ms = 10.)
# subplot(312)
# ax2 = gca()
# plot(linspace(0, rt, length(s1)), s1)
# plot(linspace(0, rt, length(s1)), s2)
# ylabel("Inputs")
# subplot(313)
# plot(linspace(0, rt, length(s)), s)
# xlabel("Time (s)")
# ylabel("Activity")
# axhline(0.5, color = "r", linestyle = "dashed")
# plot(t2*25/1000., fill(0.5, length(t2)), "k.", ms = 10.)
# plot(linspace(0, rt, length(s)), s)
# xlabel("Time (s)")
# ylabel("Activity")
# title("Normalized Activity Over Time")
# axhline(0.7, linestyle = "dashed", color="r")
# axhline(0.3, linestyle = "dashed", color="r")
# axhline(0.5, linestyle = "dashed", color="g")
# plot([i[1]/(1000/h) for i in top], fill(0.7, length(top)), "k.", ms = 10.)
# plot([i[1]/(1000/h) for i in bot], fill(0.3, length(bot)), "k.", ms = 10.)
# plot([i[1]/(1000/h) for i in nmz], fill(0.5, length(nmz)), "k.", ms = 10.)
#

# plot(t2[fw]*25/1000., fill(0.6, length(fw)), "k.", ms = 10.)
# plot(t2[fl]*25/1000., fill(0.4, length(fl)), "k.", ms = 10.)


# e_m = find(r .<= half)
# i_m = find(r .> half)
#
# te_pt = t[e_m]
# re_pt = r[e_m]
# ntd, nts = nt_diff(te_pt, re_pt, rot, quar, 150./h)
# s = ntd ./ nts
# function zscore(x)
#     m = mean(x)
#     s = std(x)
#     x1 = x .- m
#     x2 = x1 ./s
#     return x2
# end
# sd = s2 .- s1
# ss = s1 .+ s2
# sA = sd ./ ss
# #downsample until the input is the same size as the activity vector
# sA_S = downsample(sA, s)
#
# sA_S2z = zscore(sA_S)
# s2z = zscore(s)
#
# plot(sA_S2z, label = "Network Input")
# plot(s2z, label = "Network Activity")
#
# subplot(211)
# plot(te_pt, re_pt, "g.", ms = 1.)
# plot([-5000, -1000], [625, 625], lw = 5)
# plot([-5000, -1000], [1875, 1875], lw = 5)
# ylabel("Neuron Index")
# xticks([])
# subplot(212)
# xlabel("Time")
# ylabel("Current (mV/ms)")
# plot(s1/10.)
# plot(s2/10.)
# xticks([])



# JIE = [20, 30, 40, 50, 60, 70, 100, 200]
# subplot(211)
# title("Brent Network Constant Input Everywhere")
# plot(JIE, ER, "r.", ms = 10., label = "Excitatory")
# plot(JIE, IR, "b.", ms = 10., label = "Inhibitory")
# ylabel("Rates")
# legend()
# subplot(212)
# plot(JIE, ECV, "r.", ms = 10., label = "Excitatory")
# plot(JIE, ICV, "b.", ms = 10., label = "Inhibitory")
# legend()
# ylabel("CVISI")
# xlabel("JIE")
