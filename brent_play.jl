#Brent_Play
include("brent.jl")
using Distributions
srand(1234)

N = 5000
half = N/2.
quar = half/2.
runtime = 500000
rt = runtime /1000.
h = 0.1
rot = round(Int64, runtime/h)
rot2 = rot*2

stdev = .25
tau_n = 200

N_TRANS = []

# for i = [.2, .23, .26, .29, .32, .35, .38, .41, .44, .47, .5, .53, .56, .59, .62, .65, .68, .71]
#
# mu = i
mu = .2

R1 = rand(Normal(0, stdev), rot2)
R2 = rand(Normal(0, stdev), rot2)

s1 = OU_Model(R1, tau_n, 0.1)
s2 = OU_Model(R2, tau_n, 0.1)

s1 .+= mu
s2 .+= mu

s1 = s1[rot:end]
s2 = s2[rot:end]

# STD = mean([std(s1), std(s2)])
# stdev/(sqrt(2)*sqrt(tau_n))
s3 = s1 .+ s2
A1 = s1./(s3)
A2 = s2./(s3)

BW = Brent_W(N, 1/8.)
CSR = sparse_rep(BW, N)

Local_Top_Rates = []
Local_Bot_Rates = []

Top_Rates = []
Bot_Rates = []

@time t, r = Brent_Network_Euler_CSR(h, runtime, CSR, BW, N, s1, s2, 20., 20., 2.)
#for the local simulation...
# @time tl, rl = Brent_Network_Local_Euler_CSR(h, runtime, CSR, BW, N, s1, 20., 20., 2.)

# e_ml = find(rl .<= half)
# i_ml = find(rl .> half)
#
# tel_pt = tl[e_ml]
# rel_pt = rl[e_ml]

e_m = find(r .<= half)
i_m = find(r .> half)

te_pt = t[e_m]
re_pt = r[e_m]

# trl = length(find(rel_pt .> quar))
# brl = length(find(rel_pt .<= quar))
# tr = length(find(re_pt .> quar))
# br = length(find(re_pt .<= quar))
# push!(Local_Top_Rates, trl)
# push!(Local_Bot_Rates, brl)
# push!(Top_Rates, tr)
# push!(Bot_Rates, br)















ti_pt = t[i_m]
ri_pt = r[i_m]

ntd, nts = nt_diff(te_pt, re_pt, rot, quar, 150./h)
s = ntd ./ nts
flags, times = WLD_01(s, -.333, .333)

d = convert(Array{Float64}, diff(1500/10000. .* times))
cvd = cv(d)
LP = .3
dx = []
for i in d
    if i > LP
        push!(dx, i)
    end
end
dx = convert(Array{Float64}, dx)
cvdlp = cv(dx)

println("CVD is: $(cvd)")
println("CVD_LP is: $(cvdlp)")

plt[:hist](d, 50)

# push!(N_TRANS, length(flags))
# end


# s3 = s1 .+ s2 add signals together
# A1 = s1./(s3) get ratio of both signals to the sum
# A2 = s2./(s3)

ma = mean(A2) #mean of A/A+B
A2A = A2 .- ma #center A/A+B

rx = std(s)/std(A2) #rescaling factor for spikes vs inputs
A2A *= rx #rescale A/A+B
A2A += ma #put the mean back in

xtime = linspace(0, rot/10000., length(A2A))

# subplot(311)
# title("Random Network Tracks Competitive Inputs")
# plot(te_pt./10000, re_pt, "g.", ms = 1.)
# xticks([])
# yticks([625, 1875.], ["Pool 1", "Pool 2"])
# subplot(312)
# ax1 = gca()
# plot(s)
# axhline(0.5, linestyle = "dashed", color = "k")
# xticks([])
# yticks([])
# yticks([])
# ylabel("Network Activity")
# subplot(313)
# plot(xtime, A2)
# axhline(0.5, linestyle = "dashed", color = "k")
# yticks([])
# ylabel("Input Activity")
#
#
#
# println("variance of A: $(var(A2))")
# println("variance of s: $(var(s))")
# println("var(s)/var(a): $(var(s)/var(A2))")
# println("sqrt(runtime): $(sqrt(length(s)))")



function zscore(x)
    m = mean(x)
    s = std(x)
    x1 = x .- m
    x2 = x1 ./s
    return x2
end

A2z = zscore(A2)
s2z = zscore(s)
subplot(211)
plot(s2z)
subplot(212)
plot(A2z)



#difference/sum for inputs
sd = s2 .- s1
ss = s1 .+ s2
sA = sd ./ ss
#downsample until the input is the same size as the activity vector
sA_S = downsample(sA, s)

sA_S2z = zscore(sA_S)
s2z = zscore(s)
