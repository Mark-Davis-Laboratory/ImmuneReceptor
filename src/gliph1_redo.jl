module ImmuneReceptor

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

# ----------------------------------------------------------------------------------------------- #

using ProgressMeter: @showprogress

using Nucleus

# to do:
# -> add index tracking to the motif function ✅
# -> add a way to check for index / window overlap (done + updated) ✅
# -> add implementation of blosum function in the make edges
# -> combine blosum functions into one ✅
# -> create function to make like a nice summary file ✅
# -> implement fisher's exact test for cluster significance
# -> edit scoring functions to either consistently use vector of dictionary or tuple
# -> structure every scoring function to have basically the same structured output so there can be one function to read them all

# =============================================================================================== #
# Reading
# =============================================================================================== #

function is_t_gene(::Any)

    false

end

function is_t_gene(st::AbstractString)

    startswith(st, 'T')

end

function is_cdr3(st)

    st[1] == 'C' && st[end] == 'F'

end

function load_cdr3s(csvpath)
    df = CSV.read(csvpath, DataFrame)
    chain_keep_a = df.chain .∈ Ref(["TRG", "TRD", "TRA", "TRB"])
    chain_keep_b = .!(df.chain .∈ Ref(["None", "none", "NA", "na", " ", "Multi"]))
    cdr3_keep  = .!(df.cdr3 .∈ Ref(["None", "none", "NA", "na", " "]))
    c_check = startswith.(df.cdr3, "C")
    f_check = endswith.(df.cdr3, "F")

    mask = c_check .& chain_keep_a .& chain_keep_b .& cdr3_keep .& f_check

    filtered_df = df[mask, :]

    filtered_df_2 = transform(groupby(filtered_df, :cdr3), nrow => :duplicate_count)
    return filtered_df_2
end

function load_cdr3s2(csvpath)
    df = CSV.read(csvpath, DataFrame)
    cdr3_keep  = .!(df.cdr3 .∈ Ref(["None", "none", "NA", "na", " "]))
    c_check = startswith.(df.cdr3, "C")
    f_check = endswith.(df.cdr3, "F")

    mask = c_check .& cdr3_keep .& f_check

    filtered_df = df[mask, :]

    filtered_df_2 = transform(groupby(filtered_df, :cdr3), nrow => :duplicate_count)
    return filtered_df_2
end

function load_cdr3_fasta(path::AbstractString)
    cdr3s = String[]
    open(FASTA.Reader, path) do reader
        for record in reader
            push!(cdr3s, String(sequence(record)))
        end
    end
    return cdr3s
end

# =============================================================================================== #
# Writing
# =============================================================================================== #

function make_trace(an_)

    um = 1000000

    Dict(
        "name" => "Naive",
        "type" => "histogram",
        "histnorm" => "probability",
        "x" => lastindex(an_) <= um ? an_ : rand(an_, um),
    )

end

function writ(ht, st, a1_, a2_)

    Nucleus.Plotly.writ(
        ht,
        (make_trace(a1_), make_trace(a2_)),
        Dict(
            "yaxis" => Dict("title" => Dict("text" => "Probability")),
            "xaxis" => Dict("title" => Dict("text" => st)),
        ),
    )

end

# =============================================================================================== #
# MISC
# =============================================================================================== #

function make_dictionary(an_)

    Dict(an_[nd] => nd for nd in eachindex(an_))

end

# =============================================================================================== #
# CDR3
# =============================================================================================== #

function make_hamming_distance(s1, s2)

    sum(s1[nd] != s2[nd] for nd in eachindex(s1))

end

const blosum62 = Dict{Tuple{Char,Char},Int}(
    ('A','A')=>4,  ('A','R')=>-1, ('A','N')=>-2, ('A','D')=>-2, ('A','C')=>0,
    ('A','Q')=>-1, ('A','E')=>-1, ('A','G')=>0,  ('A','H')=>-2, ('A','I')=>-1,
    ('A','L')=>-1, ('A','K')=>-1, ('A','M')=>-1, ('A','F')=>-2, ('A','P')=>-1,
    ('A','S')=>1,  ('A','T')=>0,  ('A','W')=>-3, ('A','Y')=>-2, ('A','V')=>0,

    ('R','R')=>5,  ('R','N')=>0,  ('R','D')=>-2, ('R','C')=>-3, ('R','Q')=>1,
    ('R','E')=>0,  ('R','G')=>-2, ('R','H')=>0,  ('R','I')=>-3, ('R','L')=>-2,
    ('R','K')=>2,  ('R','M')=>-1, ('R','F')=>-3, ('R','P')=>-2, ('R','S')=>-1,
    ('R','T')=>-1, ('R','W')=>-3, ('R','Y')=>-2, ('R','V')=>-3,

    ('N','N')=>6,  ('N','D')=>1,  ('N','C')=>-3, ('N','Q')=>0,  ('N','E')=>0,
    ('N','G')=>0,  ('N','H')=>1,  ('N','I')=>-3, ('N','L')=>-3, ('N','K')=>0,
    ('N','M')=>-2, ('N','F')=>-3, ('N','P')=>-2, ('N','S')=>1,  ('N','T')=>0,
    ('N','W')=>-4, ('N','Y')=>-2, ('N','V')=>-3,

    ('D','D')=>6,  ('D','C')=>-3, ('D','Q')=>0,  ('D','E')=>2,  ('D','G')=>-1,
    ('D','H')=>-1, ('D','I')=>-3, ('D','L')=>-4, ('D','K')=>-1, ('D','M')=>-3,
    ('D','F')=>-3, ('D','P')=>-1, ('D','S')=>0,  ('D','T')=>-1, ('D','W')=>-4,
    ('D','Y')=>-3, ('D','V')=>-3,

    ('C','C')=>9,  ('C','Q')=>-3, ('C','E')=>-4, ('C','G')=>-3, ('C','H')=>-3,
    ('C','I')=>-1, ('C','L')=>-1, ('C','K')=>-3, ('C','M')=>-1, ('C','F')=>-2,
    ('C','P')=>-3, ('C','S')=>-1, ('C','T')=>-1, ('C','W')=>-2, ('C','Y')=>-2,
    ('C','V')=>-1,

    ('Q','Q')=>5,  ('Q','E')=>2,  ('Q','G')=>-2, ('Q','H')=>0,  ('Q','I')=>-3,
    ('Q','L')=>-2, ('Q','K')=>1,  ('Q','M')=>0,  ('Q','F')=>-3, ('Q','P')=>-1,
    ('Q','S')=>0,  ('Q','T')=>-1, ('Q','W')=>-2, ('Q','Y')=>-1, ('Q','V')=>-2,

    ('E','E')=>5,  ('E','G')=>-2, ('E','H')=>0,  ('E','I')=>-3, ('E','L')=>-3,
    ('E','K')=>1,  ('E','M')=>-2, ('E','F')=>-3, ('E','P')=>-1, ('E','S')=>0,
    ('E','T')=>-1, ('E','W')=>-3, ('E','Y')=>-2, ('E','V')=>-2,

    ('G','G')=>6,  ('G','H')=>-2, ('G','I')=>-4, ('G','L')=>-4, ('G','K')=>-2,
    ('G','M')=>-3, ('G','F')=>-3, ('G','P')=>-2, ('G','S')=>0,  ('G','T')=>-2,
    ('G','W')=>-2, ('G','Y')=>-3, ('G','V')=>-3,

    ('H','H')=>8,  ('H','I')=>-3, ('H','L')=>-3, ('H','K')=>-1, ('H','M')=>-2,
    ('H','F')=>-1, ('H','P')=>-2, ('H','S')=>-1, ('H','T')=>-2, ('H','W')=>-2,
    ('H','Y')=>2,  ('H','V')=>-3,

    ('I','I')=>4,  ('I','L')=>2,  ('I','K')=>-3, ('I','M')=>1,  ('I','F')=>0,
    ('I','P')=>-3, ('I','S')=>-2, ('I','T')=>-1, ('I','W')=>-3, ('I','Y')=>-1,
    ('I','V')=>3,

    ('L','L')=>4,  ('L','K')=>-2, ('L','M')=>2,  ('L','F')=>0,  ('L','P')=>-3,
    ('L','S')=>-2, ('L','T')=>-1, ('L','W')=>-2, ('L','Y')=>-1, ('L','V')=>1,

    ('K','K')=>5,  ('K','M')=>-1, ('K','F')=>-3, ('K','P')=>-1, ('K','S')=>0,
    ('K','T')=>-1, ('K','W')=>-3, ('K','Y')=>-2, ('K','V')=>-2,

    ('M','M')=>5,  ('M','F')=>0,  ('M','P')=>-2, ('M','S')=>-1, ('M','T')=>-1,
    ('M','W')=>-1, ('M','Y')=>-1, ('M','V')=>1,

    ('F','F')=>6,  ('F','P')=>-4, ('F','S')=>-2, ('F','T')=>-2, ('F','W')=>1,
    ('F','Y')=>3,  ('F','V')=>-1,

    ('P','P')=>7,  ('P','S')=>-1, ('P','T')=>-1, ('P','W')=>-4, ('P','Y')=>-3,
    ('P','V')=>-2,

    ('S','S')=>4,  ('S','T')=>1,  ('S','W')=>-3, ('S','Y')=>-2, ('S','V')=>-2,

    ('T','T')=>5,  ('T','W')=>-2, ('T','Y')=>-2, ('T','V')=>0,

    ('W','W')=>11, ('W','Y')=>2,  ('W','V')=>-3,

    ('Y','Y')=>7,  ('Y','V')=>-1,

    ('V','V')=>4
)

function make_blosum_score(s1, s2)

    @assert length(s1) == length(s2)

    mismatch_indices = Int[]

    for nd in eachindex(s1)
        if s1[nd] != s2[nd]
            push!(mismatch_indices, nd)
        end
    end

    if isempty(mismatch_indices) || length(mismatch_indices) > 1
            return length(mismatch_indices)
    end

    nd = mismatch_indices[1]

    a1 = s1[nd]
    a2 = s2[nd]

    scores = blosum62[(a1,a2)], blosum62[(a2,a1)]

    return min(scores)

end

function make_distance(st_, ids)

    u1 = lastindex(st_)

    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)

    po_  = Vector{Int}(undef, u2)

    i1 = 0

    @showprogress for i2 in 1:u1, i3 in (i2+1):u1

        s1 = st_[i2]
        s2 = st_[i3]

        if lastindex(s1) != lastindex(s2)

            continue

        end

        in__[i1+=1] = (Int(ids[i2]), Int(ids[i3]))
        po_[i1]     = make_hamming_distance(s1, s2)

    end

    resize!(in__, i1)
    resize!(po_,  i1)

    return in__, po_
end


# =============================================================================================== #
# Motif
# =============================================================================================== #

# count whether motif is present at least once per cdr3 for all cdr3s

function get_motif_counts(
    motif::AbstractVector{<:AbstractString},
    cdrs::AbstractVector{<:AbstractString},
)

    motif_counts = Dict{String,Int}(m => 0 for m in motif)

    cdrs = unique(cdrs)
    cores = [c[4:(end-3)] for c in cdrs if lastindex(c) >= 7]

    for cdr in cores

        for m in motif

            if occursin('.', m)

                hit = occursin(Regex(m), cdr)

            else

                hit = occursin(m, cdr)

            end

            if hit
                motif_counts[m] += 1
            end

        end

    end

    return motif_counts

end

# makes list of motifs in cdr3s, counts them, filters based on a cutoff, then provides set count
function get_motifs(st_::AbstractVector{<:AbstractString}, min::Int, max::Int, discontiguous = false)

    mo_ = Dict{String,Int}() # make dictionary

    st_ = unique(st_)
    cores = [c[4:(end-3)] for c in st_ if lastindex(c) >= 7]
    #cores = unique(cores)

    @showprogress desc = "Finding and counting motifs..." for s1 in cores

        lastindex(s1) < 7 && continue # only keep cdr3 7 or longer

        #s2 = s1[4:(end-3)] # remove first and last 3 aa
        s2 = s1
        lind = lastindex(s2)

        seen = Set{String}()

        for um in min:max

            um > lastindex(s2) && continue # make sure cdr3 is longer than the motif size

            i1 = 0
            i2 = i1 + um - 1

            for i1 in 1:(lind - um + 1)

                i2 = i1 + um - 1

                m = s2[i1:i2]

                push!(seen, m)

            end

        end

        if discontiguous == true
            for span in (4, 5)

                lastindex(s2) < span && continue

                for i in 1:(lastindex(s2) - span + 1)

                    w = s2[i:i+span-1]

                    for g in 2:(span-1)
                        p = string(w[1:g-1], ".", w[g+1:end])
                        push!(seen, p)

                    end
                end
            end
        end

        for m in seen
            mo_[m] = get!(mo_, m, 0) + 1
        end

    end

    mo2_ = Dict(m => num for (m, num) in mo_ if num >= 3)
    return mo2_

end

function find_significant_motifs(motifs, cdrs1, cdrs2, nsim, ove_cutoff)

    Random.seed!(3)

    motifs_list = collect(keys(motifs))
    L = length(motifs_list)
    counts_orig = [motifs[m] for m in motifs_list]

    counts_sim = Array{Int32}(undef, nsim, L)

    # decide whether to use regex or not
    booler = Vector{Union{Nothing,Regex}}(undef, L)

    for j in 1:L

        m = motifs_list[j]

        if occursin('.', m)
            booler[j] = Regex(m)

        else
            booler[j] = nothing

        end
    end

    significant_motifs = Dict{String,Float64}()

    cdrs2 = unique(cdrs2)
    cdrs1 = unique(cdrs1)

    cores2 = [c[4:(end-3)] for c in cdrs2 if lastindex(c) >= 7]
    cores2 = unique(cores2)

    @showprogress desc = "Calculating significant motifs..." for i in 1:nsim

        countz = zeros(Int, L)

        random_cdr_cores = sample(cores2, length(cdrs1); replace=true, ordered=false)

        for cdr in random_cdr_cores

            for j in 1:L
                    if booler[j] === nothing
                        hit = occursin(motifs_list[j], cdr)
                    else
                        hit = occursin(booler[j], cdr)
                    end

                    if hit
                        countz[j] += 1
                    end
                end

        end

        counts_sim[i, :] = countz

    end

    # everything above this line works

    for (index, m) in enumerate(motifs_list)

        if ove_cutoff == true

            obs = counts_orig[index]
            mean_sim = mean(@view counts_sim[:, index])

            # turboGLIPH-style: if expected is 0, set OVE tiny (not Inf)
            if mean_sim > 0
                ove = obs / mean_sim
            else
                ove = 1e-12
            end

            if obs < 2
                continue

            elseif (obs == 2 && ove >= 1000) ||
                   (obs == 3 && ove >= 100)  ||
                   (obs >= 4 && ove >= 10)

                wins = count(x -> x >= obs, @view counts_sim[:, index])
                p_val = (wins + 1) / (nsim + 1)

                if p_val <= 0.001
                    significant_motifs[m] = get!(significant_motifs, m, 0.0) + p_val
                end
            end
        end

        if ove_cutoff != true
            wins = count(x -> x >= counts_orig[index], counts_sim[:, index])
            p_val = (wins + 1) / (nsim + 1)
            print(p_val)
            if p_val <= 0.001
                significant_motifs[m] = get!(significant_motifs, m, 0.0) + p_val
            end
        end

    end

    println(string("Number of significant motifs identified:", length(significant_motifs)))
    return significant_motifs

end


#_________#

# takes list of cdrs and checks for a motif, making pairs
function make_motif_pairs(st_::AbstractVector{<:AbstractString}, motif::AbstractString, ids)

    u1 = lastindex(st_)
    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)
    po_  = Vector{Int}(undef, u2)

    i1 = 0

    booler = nothing

    if occursin('.', motif)
        booler = Regex(motif)
    end

    for i2 in 1:u1, i3 in (i2+1):u1

        s1, s2 = st_[i2], st_[i3]
        core1 = lastindex(s1) < 7 ? "" : s1[4:(end-3)]
        core2 = lastindex(s2) < 7 ? "" : s2[4:(end-3)]

        has1 = false
        has2 = false

        if booler === nothing
            has1 = occursin(motif, core1)
            has2 = occursin(motif, core2)

        else
            has1 = occursin(booler, core1)
            has2 = occursin(booler, core2)

        end

        if has1 && has2
            i1 += 1
            in__[i1] = (Int(ids[i2]), Int(ids[i3]))   # <-- only change
            po_[i1]  = 1
        end
    end

    resize!(in__, i1)
    resize!(po_,  i1)

    return in__, po_
end

function make_motif_pairs_new(st_::AbstractVector{<:AbstractString}, motif::AbstractString)
    u1 = lastindex(st_)
    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)
    po_  = Vector{Int}(undef, u2)

    i1 = 0

    booler = nothing

    if occursin('.', motif)
        booler = Regex(motif)
    end

    for i2 in 1:u1, i3 in (i2+1):u1

        s1, s2 = st_[i2], st_[i3]
        core1 = lastindex(s1) < 7 ? "" : s1[4:(end-3)]
        core2 = lastindex(s2) < 7 ? "" : s2[4:(end-3)]

        has1 = false
        has2 = false

        if booler === nothing

            has1 = occursin(motif, core1)
            has2 = occursin(motif, core2)

        else

            has1 = occursin(booler, core1)
            has2 = occursin(booler, core2)

        end

        if has1 && has2

            nd1 = nothing
            nd2 = nothing

            if booler === nothing

                nd1 = findfirst(motif, core1)
                nd2 = findfirst(motif, core2)

            else

                nd1 = findfirst(booler, core1)
                nd2 = findfirst(booler, core2)

            end

            a1, a2 = first(nd1), last(nd1)
            b1, b2 = first(nd2), last(nd2)

            if max(a1, b1) ≤ min(a2, b2) + 3
                i1 += 1
                in__[i1] = (i2, i3)
                po_[i1] = 1

            end
        end
    end

    resize!(in__, i1)
    resize!(po_, i1)

    return in__, po_

end


# =============================================================================================== #
# Graphing
# =============================================================================================== #

function make_vertices!(g, cdrs)
    isblank(x) = x === nothing || x === missing || x in ("None", "none", "NA", "na", "Na", " ", "  ", "")

    for row in 1:nrow(cdrs)

        label = String(cdrs.cdr3[row])

        add_vertex!(g, label, Dict(:index => row)) # add vertex

        # check for vgenes, etc. and store if they exist
        if !isblank(cdrs.vgene[row])
            g[label][:vgene] = cdrs.vgene[row]
        end

        try
            if !isblank(cdrs.jgene[row])
                g[label][:jgene] = cdrs.jgene[row]
            end
        catch
        end

        try
            if !isblank(cdrs.dgene[row])
                g[label][:dgene] = cdrs.dgene[row]
            end
        catch
        end

        if !isblank(cdrs.donor[row])
            g[label][:donor] = cdrs.donor[row]
        end

    end

    return g # return graph

end

function make_vertices_indexed!(g, cdrs)
    isblank(x) = x === nothing || x === missing || x in ("None", "none", "NA", "na", "Na", " ", "  ", "")

    if :row_index ∉ names(cdrs) # make sure rows are indexed
        cdrs.row_index = collect(1:nrow(cdrs))
    end

    for row in 1:nrow(cdrs)

        label = Int(cdrs.row_index[row])
        cdr3_aa = String(cdrs.cdr3[row])

        if !has_vertex(g, label)
            add_vertex!(g, label, Dict(:cdr3 => cdr3_aa))
        end

        # check for vgenes, etc. and store if they exist
        if !isblank(cdrs.vgene[row])
            g[label][:vgene] = cdrs.vgene[row]
        end

        try
            if !isblank(cdrs.jgene[row])
                g[label][:jgene] = cdrs.jgene[row]
            end
        catch
        end

        try
            if !isblank(cdrs.dgene[row])
                g[label][:dgene] = cdrs.dgene[row]
            end
        catch
        end

        try
            if !isblank(cdrs.donor[row])
                g[label][:donor] = cdrs.donor[row]
            end
        catch
        end

    end

    return g # return graph

end

# make a global edge w/ annotation of hamming disrance
function add_global_edge!(g, u, v, distance)
    add_edge!(g, u, v, Dict(:distance => distance))
end

# add a local edge w/ annotation of motif and the pval (from get_significant_motifs)
function add_local_edge!(g, u, v, motif::String, pval)
    add_edge!(g, u, v, Dict(
        :motifs => [motif],
        :motif_pvals => [pval],
    ))
end

# initialises the graph and makes edges
function make_edges(cdrs, motifs, isglobal, islocal)

    cdrs.cdr3 = String.(cdrs.cdr3) # ensure is string

    g = MetaGraph(
        Graph();
        label_type=Int,
        vertex_data_type=Dict{Symbol,Any},
        edge_data_type=Dict{Symbol,Any},
    )

    g = make_vertices_indexed!(g, cdrs)

    if isglobal == true

        pairs, dists = make_distance(cdrs.cdr3, cdrs.row_index)

        @showprogress desc = "Making edges..." for (index, (u, v)) in enumerate(pairs)

            d = dists[index]
            d <= 1 || continue

            if haskey(g, u, v)
                g[u, v][:distance] = d

            else
                add_global_edge!(g, u, v, d)

            end
        end

    end


    if islocal == true

       @showprogress desc = "Making edges..." for (motif, pval) in motifs

           pairs, _ = make_motif_pairs(cdrs.cdr3, motif, cdrs.row_index)

           for (u, v) in pairs

                if haskey(g, u, v)

                    if !haskey(g[u, v], :motifs)
                        g[u, v][:motifs] = String[]

                    end

                    if !haskey(g[u, v], :motif_pvals)
                        g[u, v][:motif_pvals] = Float64[]
                    end

                    push!(g[u, v][:motifs], motif)
                    push!(g[u, v][:motif_pvals], pval)

                else
                    add_local_edge!(g, u, v, motif, pval)
                end


            end
        end
    end # end of local edges

    return g

end


# =============================================================================================== #
# Clusters / Significance
#
# TO DO: multiple comparisons correction?
# =============================================================================================== #

# length score
# for each group, resample 1000 groups of group size n and measure cdr3 length distribution.

using StatsBase: mode

function find_length_pvals(g, cdrs2, sim_depth)
    clusters = connected_components(g)

    modes = Vector{Int}()
    props = Vector{Float64}()

    # get most common length for each cluster
    for cluster in clusters

        labels = label_for.(Ref(g), cluster)

        cluster_lengths = Int[]

        for lbl in labels

            haskey(g[lbl], :cdr3) || continue
            push!(cluster_lengths, length(String(g[lbl][:cdr3])))

        end

        m = mode(cluster_lengths)
        push!(modes, m)
        push!(props, count(==(m), cluster_lengths) / length(cluster_lengths))

    end

    # randomly pick n sequences from the null distribution and compare to
    cluster_sizes = [length(cluster) for cluster in clusters]
    p_vals = Float64[]

    for (index, size) in enumerate(cluster_sizes)

        sim_props = Float64[]

        for i in 1:sim_depth
            random_cdrs = sample(cdrs2, size; replace=true, ordered=false)
            cluster_lengths = [length(s) for s in random_cdrs]
            push!(sim_props, (count(==(modes[index]), cluster_lengths)) / size)
        end

        # count number of times sim beats data?
        p_val = ((count(>=(props[index]), sim_props)) + 1) / (sim_depth + 1)
        push!(p_vals, p_val)

    end

    return p_vals

end

function score_lengths(g, cdrs2) # adds length scores to the graph

    clusters = connected_components(g)
    p_vals = find_length_pvals(g, cdrs2)

    for (index, cluster) in enumerate(clusters)

        for vertex in cluster
            lbl = label_for(g, vertex)
            g[lbl][:length_pval] = p_vals[index]
        end

    end

    return p_vals

end


#_______________#

# v gene score
# hyper geomtric or parameteric when over 200 sequences in one group

using Distributions

function score_vgene(g, sim_depth)
    clusters = connected_components(g)

    # make vector w/ dictionary of vgenes
    counts = Vector{Dict{String,Int}}(undef, length(clusters))
    totals = Dict{String,Int}()
    sizes = Vector{Int}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)

        di = Dict{String,Int}()
        sizes[index] = length(cluster)

        for vertex in cluster

            lbl = label_for(g, vertex) # uses indices now
            haskey(g[lbl], :vgene) || continue
            vgene = g[lbl][:vgene]

            di[vgene] = get!(di, vgene, 0) + 1
            totals[vgene] = get!(totals, vgene, 0) + 1

        end
        counts[index] = di
    end

    # actually do p-vals after counting
    v_all = [g[cdr][:vgene] for cdr in labels(g) if haskey(g[cdr], :vgene)]
    cluster_pvals = Vector{Dict{String,Float64}}(undef, length(clusters))

    for (index, di) in enumerate(counts)

        p_vals = Dict{String,Float64}()

        if sizes[index] <= 200 # for less than 200 in a cluster

            for (vgene, count_) in di
                distribution = Hypergeometric(totals[vgene], sum(sizes) - totals[vgene], sizes[index])
                p = ccdf(distribution, count_ - 1)
                p_vals[vgene] = p
            end

        else # for when > 200 in a cluster

            sim_wins = Dict{String,Int}()

            for i in 1:sim_depth # do sim_depth random pulls from the samples and count vgenes

                # sizes[index] has cluster size - stored in earlier loop
                random_vgenes = sample(v_all, sizes[index]; replace=true, ordered=false)

                # go through every vgene
                for (vgene, orig_count) in di

                    count_ = count(x -> x == vgene, random_vgenes)

                    if count_ >= orig_count
                        sim_wins[vgene] = get!(sim_wins, vgene, 0) + 1
                    else
                        sim_wins[vgene] = get!(sim_wins, vgene, 0) + 0
                    end

                end

            end

            # go through sim wins + gen new dictionary

            for (vgene, wins) in sim_wins
                p = (wins + 1) / (sim_depth + 1)
                p_vals[vgene] = p
            end

        end

        val, key = findmin(p_vals)   # val = pval (Float64), key = vgene (String)
        cluster_pvals[index] = Dict(key => (val * length(p_vals)))

    end

    # add cluster vgene pvals to graph
    for (index, cluster) in enumerate(clusters)

        for vertex in cluster
            lbl = label_for(g, vertex)
            g[lbl][:vgene_pval] = cluster_pvals[index]
        end

    end

    return cluster_pvals

end

#_______________#

# hla score

# for each CDR, get donor
# check hla file for shared hlas
# do enrichment exactly as done for v-genes

function allele_group(allele::String)
    m = match(r"^(HLA-)?([A-Z]{1,4}\d?)", allele)
    return m === nothing ? allele : m.captures[end]
end

function score_hla(g)
    clusters = connected_components(g, hla_df)

    # make dictionaries and vectors to store stuff in ✅
    counts = Vector{Dict{String,Int}}(undef, length(clusters))

    totals = Dict{String,Int}()
    sizes = Vector{Int}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters) # count things ✅

        di = Dict{String,Int}()
        sizes[index] = length(cluster)

        for vertex in cluster

            donor = g[label_for(g, vertex)][:donor]

            row_index = findfirst(==(donor), hla_df.donor)
            row_index === nothing && continue # skip typing for donors without HLA type
            row = hla_df[row_index, 2:end]

            for allele in row

                if ismissing(allele)

                    continue

                end

                a = String(allele) # needs to be string

                totals[allele] = get(totals, allele, 0) + 1 # store total counts for distribution to compare cluster to
                di[allele] = get!(di, allele, 0) + 1

            end


        end

        counts[index] = di # store individual counts for each cluster

    end

    cluster_pvals = Vector{Dict{String,Float64}}(undef, length(clusters))

    for (index, di) in enumerate(counts)

        counts[index]
        p_vals = Dict{String,Float64}()

        for (hla_, count_) in di

            distribution = Hypergeometric(totals[index][hla_], sum(sizes) - totals[index][hla_], sizes[index])
            p = ccdf(distribution, count_ - 1)
            p_vals[hla_] = p

        end

        cluster_pvals[index] = p_vals

    end

    # generate table

    pval_df = DataFrame(cluster=Int[], allele=String[], allele_group=String[], pval=Float64[])

    for (cluster, d) in pairs(cluster_pvals)
        for (allele, p) in d
            group = allele_group(allele)
            push!(pval_df, (
                cluster = cluster,
                allele  = allele,
                group   = group,
                pval    = p
            ))
        end
    end

    return pval_df
    # heatmap_df = unstack(df, :cluster, :allele, :pval)

end


# motif scoring

using Graphs

function motif_scoring(g)

    clusters = connected_components(g)
    cluster_motifs = Vector{Vector{Tuple{String,Float64}}}(undef, length(clusters))

    for (index, cluster) in enumerate(clusters)

        in_cluster = falses(nv(g))
        in_cluster[cluster] .= true

        pairs = Tuple{String,Float64}[]

        for u in cluster

            for v in neighbors(g, u)

                if in_cluster[v] && u < v

                    lu = label_for(g, u)
                    lv = label_for(g, v)

                    if haskey(g, lu, lv) && haskey(g[lu, lv], :motifs) && haskey(g[lu, lv], :motif_pvals)

                        motifs = g[lu, lv][:motifs]
                        pvals  = g[lu, lv][:motif_pvals]

                        for ind in eachindex(motifs)
                            push!(pairs, (motifs[ind], pvals[ind]))

                        end

                    end
                end
            end
        end

        cluster_motifs[index] = pairs

    end

    return cluster_motifs

end


# =============================================================================================== #
# Summarising
# =============================================================================================== #
#

# create function to report for each cluster:
# -> number of members
# -> number of unique donors and unique clones (unique alpha and beta?) - less of a priority
# -> length of cdr3 p-val
# -> motif p-val
# -> v-gene pval
# -> hla p-vals
# -> list of probable hla types

# length -> vector of p-vals
# indexed dictionary with correct p-val of most significant v-gene

function summarize_data(g, vgene_pvals, length_pvals, cluster_motifs)

    # get cluster sizes

    clusters = connected_components(g)

    # things to count:
    # unique members, unique cdr3s, (later - unique clones)
    # find most significant motif p-val

    # make vector w/ dictionary of vgenes
    cluster_sizes = Vector{Dict{String,Int}}(undef, length(clusters))
    totals = Dict{String,Int}()
    cluster_sizes = Vector{Int}(undef, length(clusters))
    donor_sizes = Vector{Int}(undef, length(clusters))
    cdr_sizes = Vector{Int}(undef, length(clusters))
    length_pval = Vector{Float64}(undef, length(clusters))

    # get cluster sizes and no. donors
    for (index, cluster) in enumerate(clusters)

        cluster_sizes[index] = length(cluster)

        donors = String[]
        cdrs = String[]

        for vertex in cluster

            lbl = label_for(g, vertex)

            cdr3 = g[lbl][:cdr3]
            push!(cdrs, String(cdr3))

            if haskey(g[lbl], :donor)

                d = g[lbl][:donor]
                d === missing && continue
                push!(donors, String(d))

            end

        end

        donor_sizes[index] = length(unique(donors))
        cdr_sizes[index] = length(unique(cdrs))
    end

    # get v-genes
    vgenes = [first(keys(i)) for i in vgene_pvals]
    v_pvals  = [first(values(i)) for i in vgene_pvals]

    # get motifs
    top_motifs = Tuple{String,Float64}[]

    for pairs in cluster_motifs
        if isempty(pairs)

            push!(top_motifs, ("", NaN))

        else

            pvals = last.(pairs)
            ind = argmin(pvals)
            push!(top_motifs, pairs[ind])

        end

    end

    motifs = first.(top_motifs)
    m_pvals  = last.(top_motifs)

    df = DataFrame(
        group = collect(1:length(clusters)),
        no_member = cluster_sizes,
        no_cdrs = cdr_sizes,
        no_donor = donor_sizes,
        length_pval = length_pvals,
        vgene_pval = v_pvals,
        vgene = vgenes,
        motif_pval = m_pvals,
        motif = motifs
    )

    return df

end

using Graphs
using DataFrames

function extract_clusters_table(g)
    clusters = connected_components(g)

    rows = NamedTuple[]

    for (group_id, cluster) in enumerate(clusters)
        for v in cluster
            lbl = label_for(g, v)  # row_index label (your vertex label)

            cdr3  = haskey(g[lbl], :cdr3)  ? String(g[lbl][:cdr3])  : ""
            vgene = haskey(g[lbl], :vgene) ? String(g[lbl][:vgene]) : ""
            jgene = haskey(g[lbl], :jgene) ? String(g[lbl][:jgene]) : ""
            dgene = haskey(g[lbl], :dgene) ? String(g[lbl][:dgene]) : ""
            donor = haskey(g[lbl], :donor) ? String(g[lbl][:donor]) : ""

            push!(rows, (
                group  = group_id,
                vertex_id = v,
                row_index  = lbl,
                cdr3   = cdr3,
                donor  = donor,
                vgene  = vgene,
                jgene  = jgene,
                dgene  = dgene,
            ))
        end
    end

    members = DataFrame(rows)

    summary = combine(groupby(members, :group),
        nrow => :no_member,
        :cdr3  => (x -> length(unique(x))) => :no_unique_cdr3,
        :donor => (x -> length(unique(filter(!isempty, x)))) => :no_donor,
        :vgene => (x -> length(unique(filter(!isempty, x)))) => :no_vgene,
    )

    return (members = members, summary = summary)
end
