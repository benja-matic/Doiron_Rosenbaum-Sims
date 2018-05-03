function emptiness(x, funk, error_code) #checks if an empty set got populated
  if length(x) > 0
    return funk(x)
  else
    return error_code
  end
end

function cv(isi)
  SD = std(isi)
  return SD/mean(isi)
end

function Brent_W(N, k, Jee, Jie, Jei, Jii)

  half = round(Int64, N/2)
  sk = sqrt(k)
  W = zeros(N, N)

  Jee /= sk
  Jie /= sk
  Jei /= sk
  Jii /= sk

  #Excitatory Neurons
  for i = 1:half
    ee_inds = rand(1:half, k)
    ei_inds = rand(half+1:N, k)
    for j in eachindex(ee_inds)
      W[i, ee_inds[j]] += Jee
      W[i, ei_inds[j]] += Jei
    end
  end

  #Inhibitory Neurons
  for i = half+1:N
    ie_inds = rand(1:half, k)
    ii_inds = rand(half+1:N, k)
    for j in eachindex(ie_inds)
      W[i, ie_inds[j]] += Jie
      W[i, ii_inds[j]] += Jii
    end
  end

  return W
end

#indices for parse column representation
function sparse_rep(W,N)
    flat = []
    for i = 1:N
        mi = find(W[:,i])
        push!(flat, mi)
    end
    return flat
end


  #20000 thousand neurons, 10000 E and 10000 I
  #each presynaptic neuron connects to 2500 random E and 2500 random I neurons
  #sampling includes replacement; can make multiple connections
  #ee = 12.5, ie = 20., ii=ei = -50, scaled by sqrt(N)
  #cov(s(t), s(t+tau)) = exp(-tau^2/2(tau^2))
  #F(t) = sqrt(N)*bias + sigma*s(t)

function OU_Model(R, tau, h)
  smooth = zeros(length(R))
  x = rand()
  scale = sqrt(h)/tau
  leak = -h/tau
  R *= scale
  for i = 1:length(R)
    x += leak*x + R[i]
    smooth[i] = x
  end
  return smooth
end

function interpolate_spike(v2, v1, vth)
  x = (v1-v2) #slope by linear interpolation (dv/dt) = change in voltage for a single time step
  t = (vth - v2)/x #time since spike to now
  return t
end

function Brent_Network_Euler_CSR_A(h, total, CSR, W, N, s1, s2, vth, tau_m, tau_s, tau_a, g_a)

  ntotal = round(Int,total/h)
  half = round(Int64, N/2)
  quar = round(Int64, N/4)

  half_ = half + 1
  quar_ = quar + 1

  V = (rand(N)*vth) .- 60
  V_buff = V

  syn = zeros(N)
  drive = zeros(N)
  A = zeros(N)

  # drive[1:quar] = s1 *h
  # drive[quar+1:half] = s2 *h
  # drive[half+1:half+quar] = (s1 - .1)*h
  # drive[half+quar:N] = (s2 - .1)*h
  g_h = g_a * h

  p1 = [1:quar] #pool 1 E
  p2 = [quar_:half] #pool 2 E
  p3 = [half_:half+quar] #pool 1 I
  p4 = [half_+quar:N] #pool 2 I
  s3 = s1 .- (h*.005)
  s4 = s1 .- (h*.005)

  #raster
  time = Float64[0]
  raster = Float64[0]

  #Calculate drive and leak constants ahead of time
  m_leak = h/tau_m
  s_leak = h/tau_s
  a_leak = h/tau_a

  for iter = 1:ntotal

    V .+= (h*syn/tau_s) .- (V*m_leak) .- (g_h*A)

    V[p1] += s1[iter]
    V[p2] += s2[iter]
    V[p3] += s3[iter]
    V[p4] += s4[iter]

    syn .-= (syn*s_leak)
    A .-= (A*a_leak)

    #check for spikes
    vs = (V.>vth)
    vsm = sum(vs)

    #update excitatory synapses and adaptation
    if vsm > 0
      spe = find(vs)
      for j = 1:vsm
        js = spe[j]
        #interpolate spike time
        delta_h = interpolate_spike(V[js], V_buff[js], vth) #time since the spike (units of h)
        lx = exp(-delta_h*h/tau_s)
        syn[CSR[js]] += W[CSR[js], js].*lx #modify amplitude of synapse by decay since estimated time of spike
        push!(raster,js)
  	    push!(time,iter-delta_h)
        A[js] += 1.
      end
    end

    #reset after spike
    # V .-= vth*vs
    V[vs] = -65.
    V_buff = V
  end

  return time, raster
end

function nt_diff(t, r, ntotal, half, netd_binsize)

  netd_bins = collect(1:netd_binsize:ntotal)
  ntd = zeros(length(netd_bins)-1)
  nts = zeros(length(netd_bins)-1)

  for j = 2:length(netd_bins)
    tf = find(netd_bins[j-1] .<= t .< netd_bins[j])
    T = sum(r[tf] .> half)
    B = sum(r[tf] .<= half)
    NTD = T - B #################Differences in Spikes
    ntd[j-1] = NTD
    nts[j-1] = length(tf)
    #println(T+B==length(tf))
  end

  return ntd, nts
end

function shmitt_single(s)

  #Handle WTA and Normalization Exceptions
  if maximum(s) < .5
    return return ["lose", "end"], [1, length(s)]
  elseif minimum(s) >= .5
    return ["win", "end"], [1, length(s)]
  elseif (.3 <= maximum(s) <= .7) & (.3 <= minimum(s) <= .7)
    return ["draw", "end"], [1, length(s)]
  end

  times = []
  flags = []

  w = findfirst(s .> .5)
  l = findfirst(s .< .5)
  f = minimum([w, l])
  if f == w
    flag = true
  else
    flag = false
  end

  push!(times, f-1)
  push!(flags, flag)

  for i = f:length(s)
    f1 = (s[i] > .5)
    if f1 != flag
      push!(times, (i+i-1)/2.)
      push!(flags, flag)
      flag = !flag
    end
  end

  return flags, times
end

function WLD_01(s, tl, th)
  if maximum(s) < tl
    return ["lose", "end"], [1, length(s)]
  elseif minimum(s) >= th
    return ["win", "end"], [1, length(s)]
  elseif (tl <= maximum(s) <= th) & (tl <= minimum(s) <= th)
    return ["draw", "end"], [1, length(s)]
  end
  times = []
  flags = []
  if s[1] >= th
    flag = "win"
  elseif tl <= s[1] <= th
    flag = "draw"
  elseif s[1] < tl
    flag = "lose"
  end

  push!(times, 1)
  push!(flags, flag)

  s2 = s[2:end]
  for i in eachindex(s2)
    f1 = comp_01(s2[i], tl, th)
    if f1 != flag
      flag = f1
      push!(times, i+1)
      push!(flags, flag)
    end
  end
  push!(times, length(s))
  push!(flags, "end")
  return flags, times
end

function comp_01(x, tl, th)
  if x > th
    return "win"
  elseif tl <= x <= th
    return "draw"
  elseif x < tl
    return "lose"
  else
    return "weird"
  end
end

function splice_reversions(flags, times)
  w0 = findfirst(flags, "win")
  l0 = findfirst(flags, "lose")
  #assume rivalry, and that any reversions are short and fail to persist
  empezar = min(w0, l0)
  nf = [flags[empezar]]
  a = ["win", "lose"]
  t = [empezar]
  if start == w0
    f = 2
  else
    f = 1
  end
  flag = a[f]
  for i = empezar:length(flags)
    if flags[i] == a[f]
      push!(t, times[i])
      push!(nf, flags[i])
      a = circshift(a, 1)
    end
  end
  push!(t, times[end])
  return t, nf
end

function splice_flags(flags, times)
  l = find(flags.=="lose")
  w = find(flags.=="win")
  d = find(flags.=="draw")
  bot = [[((times[i]-1)*250)+1, ((times[i+1]-1)*250)] for i in l]
  top = [[((times[i]-1)*250)+1, ((times[i+1]-1)*250)] for i in w]
  nmz = [[((times[i]-1)*250)+1, ((times[i+1]-1)*250)] for i in d]
  # bot = [[times[i], times[i+1]] for i in l]*250
  # top = [[times[i], times[i+1]] for i in w]*250
  # nmz = [[times[i], times[i+1]] for i in d]*250
  bdom = emptiness([(i[2]-i[1]) for i in bot], sum, 0)
  tdom = emptiness([(i[2]-i[1]) for i in top], sum, 0)
  tnmz = emptiness([(i[2]-i[1]) for i in nmz], sum, 0)
  # bdom = emptiness([250*(i[2]-i[1]) for i in bot], sum, 0)
  # tdom = emptiness([250*(i[2]-i[1]) for i in top], sum, 0)
  # tnmz = emptiness([250*(i[2]-i[1]) for i in nmz], sum, 0)
  return top, tdom, bot, bdom, nmz, tnmz
end

#downsample A so it is the same size as B
function downsample(A, B)
  LB = length(B)
  f = Int64(floor(length(A)/LB))
  L = Int64(f*LB)
  A_ = A[1:Int64(L)] #Truncated version of A so B fits into it an even number of times
  println(length(A_))
  A_S = zeros(LB) #subsample of A
  for i = 1:LB-1
    println(i)
    A_S[i] = A_[((i-1)*f)+1]
  end
  return A_S
end

function zscore(x)
    m = mean(x)
    s = std(x)
    x1 = x .- m
    x2 = x1 ./s
    return x2
end
