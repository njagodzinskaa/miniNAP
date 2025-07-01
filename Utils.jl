function loadData(filename)

    data = h5open(filename, "r") |> fid -> read(fid, "Data")
    traces = Float64.(data["Recording_0"]["AnalogStream"]["Stream_0"]["ChannelData"])
    info = data["Recording_0"]["AnalogStream"]["Stream_0"]["InfoChannel"]

    return traces, info
end

mutable struct InfoStruct
    Label::String
end

function loadDataMat(filename)

    vars = matread(filename)
    traces = vars["dat"]

    els = split(vars["header"], "\n")
    m = match.(r"\d{2}", split(els[7], ";"))
    m = [m[x].match for x in 1:length(m)]


    info = DataFrame(Label = Vector{Any}(undef, length(m)))
    for x in 1:length(m)
        info[x, :Label] = m[x]
    end

    return traces, info
end

#---------------------------------------------------------------------------

function filterTrace(traces, fs, lowpass, highpass)

    # filter_order = 3
    # responsetype = DSP.Bandpass(lowpass, highpass; fs)
    # designmethod = Butterworth(filter_order)

        
    Threads.@threads for channel in axes(traces)[2]
        traces[:, channel] = filtfilt(digitalfilter(Bandpass(0.03,0.3),Butterworth(3)), traces[:, channel])
        # traces[:, channel] .= filt(digitalfilter(responsetype, designmethod), traces[:, channel])
    end

    return traces

end

#---------------------------------------------------------------------------

function alignSpikes(trace, spk_idx, win = 20)

    waveforms = zeros(win * 2 + 1, length(spk_idx))
    new_spk_idx = zeros(length(spk_idx))

    for (i, spk) ∈ enumerate(spk_idx)
        if (spk - win < 1) || (spk + win > length(trace))  # skip spikes near the edges
            continue
        end

        waveform = trace[spk-win:spk+win]
        neg_peak_relative = findall(waveform .== minimum(waveform))[1]
        neg_peak_absolute = spk - win + neg_peak_relative - 1

        if (neg_peak_absolute - win >= 1) && (neg_peak_absolute + win <= length(trace))
            new_spk_idx[i] = neg_peak_absolute
            waveforms[:, i] .= trace[neg_peak_absolute-win:neg_peak_absolute+win]
        end
    end

    return sort!(new_spk_idx), waveforms
end

#---------------------------------------------------------------------------

function detectSpikes(traces::Array{Float64, 2}, threshold_multiplier::Float64)

    num_traces = size(traces, 2)
    spikes = Vector{Vector{Float64}}(undef, num_traces)
    waveforms = Vector{Matrix{Float64}}(undef, num_traces)

    # Threads.@threads for t in 1:num_traces
    for t in 1:num_traces
        trace = @view traces[:, t]
        mad = median(abs.(trace .- mean(trace))) / 0.6745
        threshold = median(trace) - (threshold_multiplier * mad)
        peaks = Vector{Int}(undef, 0)
        peaks = findall(trace .< threshold)
        spikes[t], waveforms[t] = alignSpikes(trace, peaks)
    end

    return spikes, waveforms
end

#---------------------------------------------------------------------------

function eye(n)
    return Int.(diagm(ones(n)))
end

#---------------------------------------------------------------------------

function run_P(N1::Int64, N2::Int64, dt::Float64, spike_times_1::Vector{Float64}, spike_times_2::Vector{Float64})
    Nab = 0
    j = 1
    for i ∈ 1:N1
        while j <= N2
            if abs(spike_times_1[i] - spike_times_2[j]) <= dt
                Nab += 1
                break
            elseif spike_times_2[j] > spike_times_1[i]
                break
            else
                j += 1
            end
        end
    end
    return Nab
end

#---------------------------------------------------------------------------

function run_T(N1::Int64, dt::Float64, start::Float64, endv::Float64, spike_times_1::Vector{Float64})
    time_A = 2 * N1 * dt
    if N1 == 1
        if (spike_times_1[1] - start) < dt
            time_A = time_A - start + spike_times_1[1] - dt
        elseif (spike_times_1[1] + dt) > endv
            time_A = time_A - spike_times_1[1] - dt + endv
        end
    else
        for i ∈ 1:(N1-1)
            diff = spike_times_1[i+1] - spike_times_1[i]
            if diff < 2 * dt
                time_A = time_A - 2 * dt + diff
            end
        end
        if (spike_times_1[1] - start) < dt
            time_A = time_A - start + spike_times_1[1] - dt
        end
        if (endv - spike_times_1[N1]) < dt
            time_A = time_A - spike_times_1[N1] - dt + endv
        end
    end
    return time_A
end

#---------------------------------------------------------------------------

function run_sttc(N1::Int64, N2::Int64, dt::Float64, Time::Vector{Float64}, spike_times_1::Vector{Float64}, spike_times_2::Vector{Float64})
    if N1 == 0 || N2 == 0
        index = NaN
    else
        T = Time[2] - Time[1]
        TA = run_T(N1, dt, Time[1], Time[2], spike_times_1) / T
        TB = run_T(N2, dt, Time[1], Time[2], spike_times_2) / T
        PA = run_P(N1, N2, dt, spike_times_1, spike_times_2) / N1
        PB = run_P(N2, N1, dt, spike_times_2, spike_times_1) / N2
        index = 0.5 * (PA - TB) / (1 - TB * PA) + 0.5 * (PB - TA) / (1 - TA * PB)
    end
    return index
end

#---------------------------------------------------------------------------

function getSTTC(spike_times, dt, duration_s)

    Time = vec([0.0 duration_s])
    num_channels = length(spike_times)
    sttc_values = Array{Float64}(undef, binomial(num_channels, 2))

    for (idx, pair) in enumerate(combinations(1:num_channels, 2))

        # pairs = collect(combinations(1:num_channels, 2))
        # Threads.@threads for idx = 1:length(pairs)
        #     pair = pairs[idx]

        spike_times_1 = spike_times[pair[1]]
        spike_times_2 = spike_times[pair[2]]
        N1 = length(spike_times_1)
        N2 = length(spike_times_2)
        sttc_values[idx] = run_sttc(N1, N2, dt, Time, spike_times_1, spike_times_2)
    end

    adjM = getAdjM(sttc_values, num_channels)
    adjM[isnan.(adjM)] .= 0
    adjM[adjM.<0] .= abs.(adjM[adjM.<0])

    return adjM
end

#---------------------------------------------------------------------------

function getAdjM(vec::Vector, dim::Int)
    # Create a zero matrix of appropriate dimension
    mat = zeros(dim, dim)

    # Index for the vector
    idx = 1

    # Fill in the lower triangle of the matrix
    for i ∈ 2:dim
        for j ∈ 1:(i-1)
            mat[i, j] = vec[idx]
            idx += 1
        end
    end

    # The matrix is symmetric, so the upper triangle is the transpose of the lower triangle
    mat = mat + transpose(mat)

    return mat
end

#---------------------------------------------------------------------------

function probThr(spike_times, rep_num, duration_s, dt, tail)

    num_nodes = length(spike_times)
    A_new = getSTTC(spike_times, dt, duration_s)
    adjMi = get_adjMi(spike_times, dt, duration_s, rep_num)

    for i in 1:num_nodes
        for j in 1:num_nodes
            topval = quantile(adjMi[i, j, :], 1 - tail)
            if A_new[i, j] < topval
                A_new[i, j] = 0
            end
        end
    end

    return A_new
end

#---------------------------------------------------------------------------

function get_adjMi(spike_times, dt, duration_s, rep_num)
    num_nodes = length(spike_times)
    adjMi = Array{Float64, 3}(undef, num_nodes, num_nodes, rep_num)
    synth_spk = [Float64[] for _ in 1:num_nodes]
    for i ∈ 1:rep_num
        for n ∈ 1:num_nodes
            t = rand(1:duration_s)
            synth_spk[n] = sort(circshift_spikes(spike_times[n], t, duration_s))
        end
        adjMs = getSTTC(synth_spk, dt, duration_s)
        inds = LinearIndices(adjMs)
        adjMs[inds[diagind(adjMs)]] .= 0.0
        adjMi[:, :, i] .= adjMs
    end
    return adjMi
end

#---------------------------------------------------------------------------

function circular_shift(spike_times::Vector{Float64}, shift_val::Float64)
    n = maximum(spike_times)
    shift_val = mod(shift_val, n)
    return (mod.(spike_times .+ shift_val, n))
end

#---------------------------------------------------------------------------

function prune(adjMorg, pct = 0.1)
    adjM = adjMorg
    numel = length(adjM[:])
    last_idx = Int(ceil(numel * pct))
    maxA = sort(adjM[:], rev = true)[last_idx] # Get only top 10% values 
    adjM[adjM.<maxA] .= 0
    return adjM
end

#---------------------------------------------------------------------------

function ave_ctrb_lyap(A)
    B = Diagonal(ones(60))
    A[Bool.(B)] .= 1 # set diagonal to 1
    D = svd(A)
    A = A ./ (1 + maximum(D.S))
    A -= eye(60) # normalization

    W = lyapc(A, B) # Solve Lyapunov equation
    ctrb = W[Bool.(B)] # Get trace of ctrb Gramian

    return ctrb
end

#---------------------------------------------------------------------------

function circshift_spikes(spike_times, t, duration_s)
    spk_vec = spike_times .+ t
    overhang = spk_vec .> duration_s
    spk_vec[overhang] = spk_vec[overhang] .- duration_s
    return sort(spk_vec)

end

#---------------------------------------------------------------------------

function plotHeatmap(values_original, ref_el_id, stim_elec_id, annotate_flag, marker_size, cbar_title, min_val, max_val, type)

    # type = "diff"
    values = [values_original; max_val]

    hM = zeros(8, 8) .+ min_val

    # cgrads = cgrad([:teal, :white, :orange], 256)
    cgrads = cgrad([ColorSchemes.viridis[1], :white, ColorSchemes.viridis[end]], 256)

    h = plot()

    for (ctr, channel) in enumerate(chan_labels)

        if type == "diff"
            rescaled_value = values[ctr] ./ maximum([abs(min_val), abs(max_val)])
            color_value = get(cgrads, rescaled_value * 0.5 + 0.5)
            symmetric_max_val = maximum([abs(max_val) abs(min_val)])
            hM[1, 1] = symmetric_max_val
            hM[8, 1] = -symmetric_max_val
        else
            colormap = ColorSchemes.plasma
			
            # Handle NaN values: assign a default color index (e.g., 1) for NaNs
            color_id_vec = fill(1, length(values))
            not_nan_idx = findall(.!isnan.(values))
            color_id_vec[not_nan_idx] .= Int64.(round.(rescale(values[not_nan_idx], 1, 256)))

            if isnan(values[ctr])
                color_value = :gray  # or any arbitrary color for NaN
            else
                color_value = colormap[color_id_vec[ctr]]
            end
            cgrads = :plasma
            hM[1, 1] = max_val
            hM[8, 1] = min_val
        end

        y = Int64(mod(channel, 10))
        x = Int64((channel - y) / 10)
        hM[x, y] = values[ctr]

        if ctr in ref_el_id
            annotate!(x, y, "")
            continue
        end

        h = plot!([x], [y],
            seriestype = :scatter,
            markersize = marker_size,
            # markerstrokecolor=:transparent,
            markerstrokewidth = 0,
            aspect_ratio = :equal,
            color = color_value,
            # markershape=:square,
            markershape = :circle,
            xlim = (0, 9),
            ylim = (0, 9),
            axis = false,
            grid = false,
            label = "",
            legend = :outertopleft,
        )

        if ctr == stim_elec_id
            annotate!(x, y, "#")
            # h = plot!([x], [y],
            # seriestype=:scatter, markerstrokecolor=:red, markerstrokewidth=5, markersize=20, markerfillalpha=0, label="", markershape=:square)
            continue
        end

        # ----------- ANNOTATE CHANNEL NUMBER -----------
        # (This line is new:)
        annotate!(x, y, text(string(channel), :black, 8, :center))
        # ^^^ Change :black and 8 as you like
        # ----------------------------------------------
    end

    # Annotate corners and/or reference electrodes
    if annotate_flag
        annotate!(1, 1, "11")
        annotate!(1, 8, "18")
        annotate!(8, 8, "88")
        annotate!(8, 1, "81")
    end

    colorbar_labels = (min_val, max_val, 5)

    yflip!(true)
    heatmap!(hM, alpha = 0.0, color = cgrads, colorbar = true, colorbar_title = cbar_title, xlim = (0, 9), ylim = (0, 9))
    Plots.gr_cbar_width[] = 0.01
    return h
end

#---------------------------------------------------------------------------

function rescale(vector, new_min, new_max)
    old_min = minimum(vector)
    old_max = maximum(vector)
    return ((vector .- old_min) ./ (old_max - old_min)) .* (new_max - new_min) .+ new_min
end

#---------------------------------------------------------------------------

# Plot spike markers
function plotSpikeMarkers(filtered_traces, spikes, channel, threshold_multiplier)

    trace = filtered_traces[:, channel]
    mad = median(abs.(trace .- mean(trace))) / 0.6745
    hei = median(trace) - (threshold_multiplier * mad)
    p = plot(trace)

    plot!(spikes[channel] .* fs,
        zeros(length(spikes[channel])) .+ hei,
        seriestype = :scatter, xlim = (1, 100_000),
        markershape = :utriangle,
        grid = false,
        background_color_inside = :transparent,
        background = :transparent,
        legend = false,
        size = (1000, 300),
        dpi = 300,
    )

    hline!([hei])
    return p
end

#---------------------------------------------------------------------------

function ave_ctrb(A)
    A = A ./ (1 + svds(A, nsv = 1)[1].S[1])  # Matrix normalization
    T, U = schur(A)                        # Schur stability
    midMat = permutedims(U .^ 2)

    v = diag(T)
    P = repeat(diag(1 .- v * v'), 1, size(A, 1))
    values = sum(midMat ./ P, dims = 1)'
    return values
end

#---------------------------------------------------------------------------

function main(filename)
    # filename = processed_files[1]

    # traces, info = loadData(filename)
    traces, info = loadDataMat(filename)

    # channel_id = [info[i].ChannelID + 1 for i in eachindex(info)]
    # channel_id = [info[i].ChannelID + 1 for i in length(info)]

    channel_id = collect(1:60)

    # global chan_labels = [info[i].Label == "Ref" ? Int32(15) : parse(Int32, info[i].Label) for i in eachindex(info)]
    # global chan_labels = [info[i].Label == "Ref" ? Int32(15) : parse(Int32, info[i].Label) for i in 1:nrow(info)]
    global chan_labels = [label == "Ref" ? Int32(15) : parse(Int32, label) for label in info[!, :Label]]

    grounded = [ch for ch in channels if !(ch in chan_labels)]

    # Settings
    global fs = 25_000
    dt = 0.05 # time delay in seconds for STTC
    lowpass = 8_000
    highpass = 600
    threshold_multiplier = 4.0 # times median absolute deviation of the signal
    duration_s = size(traces, 1) / fs

    # Filtering (not sure if needs to be done)
    filtered_traces = filterTrace(traces, fs, highpass, lowpass)
    # filtered_traces = traces

    # Spike detection
    spikes_org, waveforms = detectSpikes(filtered_traces, threshold_multiplier) ./ fs

    # plot(waveforms[55], color = :white, linewidth = 0.1, alpha = 0.2, legend = false, grid = false)
    # plot!(mean(waveforms[55], dims = 2), color = :black, linewidth = 2)

    # Get spiking frequencies
    spk_freq = zeros(length(spikes_org))
    for i ∈ eachindex(spikes_org)
        spk_freq[i] = length(spikes_org[i]) / duration_s
    end

    id = findall(channels .== 67)[1]
    # h = plotSpikeMarkers(filtered_traces, spikes_org, id, threshold_multiplier)

    num_active_channels = sum(spk_freq .> 0.0)
    active_channel_ids = findall(spk_freq .> 0.0)

    # plotSpikeMarkers(spikes, 47, threshold_multiplier)

    spikes = spikes_org[active_channel_ids]
    A = getSTTC(spikes, dt, duration_s)
    h1 = heatmap(A, title = "Raw", aspect_ratio = :equal, axis = false, grid = false, background_color_inside = "transparent")

    rep_num = 300
    tail = 0.05
    adjMci = probThr(spikes, rep_num, duration_s, dt, tail)
    h3 = heatmap(adjMci, title = "Thresholded")

    plot(h1, h3, layout = (1, 2), aspect_ratio = :equal, grid = false, axis = false, background_color_inside = :transparent, colorbar = false)
    yflip!(colorbar = true)

    adjMci = prune(adjMci)
    heatmap(adjMci)
    yflip!()

    ctrb = ave_ctrb(adjMci)
	print("Controllability: ", ctrb)
    

    marker_size = 20

    average_ctrb = fill(minimum(ctrb), length(spikes_org))
    average_ctrb[active_channel_ids] .= ctrb

    nonactive_ids = collect(1:length(spikes_org))
    nonactive_ids = setdiff(nonactive_ids, active_channel_ids)

    min_AC = 0.0
    max_AC = maximum(ctrb)

    h1 = plotHeatmap(average_ctrb, [nonactive_ids; grounded], [], true, marker_size, "", min_AC, max_AC, "")
    title!(h1, "Controllability")
    h2 = plotHeatmap(spk_freq, nonactive_ids, [], true, marker_size, "", min_AC, max_AC, "")
    title!(h2, "Spiking frequency")
    h3 = plot(h1, h2, layout = (1, 2), markersize = 10, background_color_inside = :black, background = :black, dpi = 300)

    savefig(h3, "heatmaps.png")

    highest = rescale(average_ctrb, 0, 10) .* rescale(spk_freq, 0, 0.1)
    plotHeatmap(highest, nonactive_ids, [], false, marker_size, "", 0, 1, "")

    lowest = rescale.((1 .- average_ctrb), 0, 1) .* rescale.(spk_freq, 0, 1)
    plotHeatmap(lowest, nonactive_ids, [], false, marker_size, "", 0, 1, "")

    lowest = rescale(1 .- average_ctrb, 0, 1) .* rescale(spk_freq, 0, 1)
    low_AC_coord = channels[findall(lowest .== maximum(lowest))]

    highest = rescale(average_ctrb, 0, 1) .* rescale(spk_freq, 0, 1)
    high_AC_coord = channels[findall(highest .== maximum(highest))]

    df = DataFrame(ctrb = average_ctrb, spk_freq = spk_freq, channel = chan_labels)
    sort!(df, :ctrb, rev = true)
    CSV.write("average_controllability_2C_250701.csv", df)

    # top5 = df.ctrb[1:5] .* df.spk_freq[1:5]
    top5 = df.spk_freq[1:5]
    top_idx = minimum(findall(top5 .== maximum(top5)))

    top_coord = df.channel[top_idx]
    top_ac = df.ctrb[top_idx]
    top_fr = df.spk_freq[top_idx]

    bot_fr = df.spk_freq[df.ctrb.>1.0]
    bot5 = bot_fr[end-4:end]
    bot_idx = maximum(findall(bot_fr .== maximum(bot5)))
    bot_coord = df.channel[bot_idx]
    bot_ac = df.ctrb[bot_idx]
    bot_fr = df.spk_freq[bot_idx]

    df_row = DataFrame(recording = filename, high_coord = top_coord[1], high_AC = top_ac[1], high_fr = top_fr[1], low_coord = bot_coord[1], low_AC = bot_ac[1], low_fr = bot_fr[1])

    append!(res, df_row)

    return df
end

function cosine_similarity(A, B)
    flat_A = vec(A)
    flat_B = vec(B)
    return dot(flat_A, flat_B) / (norm(flat_A) * norm(flat_B))
end
