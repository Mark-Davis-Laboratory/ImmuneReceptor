module ImmuneReceptor

const PK = pkgdir(ImmuneReceptor)

const IN = joinpath(PK, "in")

const OU = joinpath(PK, "ou")

# ----------------------------------------------------------------------------------------------- #

using ProgressMeter: @showprogress

using Nucleus

# TO DO (in progress) - allow for analysis of multiple chains per cell, paired graphing, etc.


# =============================================================================================== #
# Reading / Misc
# =============================================================================================== #

# ✅ checked
function is_t_gene(::Any)

    false

end

# ✅ checked
function is_t_gene(st::AbstractString)

    startswith(st, 'T')

end

# ✅ checked
function is_cdr3(st)

    st[1] == 'C' && st[end] == 'F'

end

# ✅ checked
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
# VJ
# =============================================================================================== #

function make_dictionary(an_)

    Dict(an_[nd] => nd for nd in eachindex(an_))

end

function make_vj(s1_, s2_)

    u1_ = unique(s1_)

    u2_ = unique(s2_)

    U = zeros(Int, lastindex(u1_), lastindex(u2_))

    d1 = make_dictionary(u1_)

    d2 = make_dictionary(u2_)

    for nd in eachindex(s1_)

        U[d1[s1_[nd]], d2[s2_[nd]]] += 1

    end

    u1_, u2_, U

end

# =============================================================================================== #
# CDR3
# =============================================================================================== #

# ✅ checked
function make_hamming_distance(s1, s2)

    sum(s1[nd] != s2[nd] for nd in eachindex(s1))

end

# ✅ checked
function make_distance(st_)

    u1 = lastindex(st_)

    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)

    po_ = Vector{Int}(undef, u2)

    i1 = 0

    @showprogress for i2 in 1:u1, i3 in (i2+1):u1

        s1 = st_[i2]

        s2 = st_[i3]

        if lastindex(s1) != lastindex(s2)

            continue

        end

        in__[i1+=1] = i2, i3

        po_[i1] = make_hamming_distance(s1, s2)

    end

    resize!(in__, i1), resize!(po_, i1)

end

# =============================================================================================== #
# Motif
# =============================================================================================== #

# count whether motif is present at least once per cdr3 for all cdr3s

# ✅ checked
function get_motif_counts(
    motif::AbstractVector{<:AbstractString},
    cdrs::AbstractVector{<:AbstractString},
)

    motif_counts = Dict{String,Int}()

    for m in motif

        motif_counts[m] = count(cdr -> occursin(m, cdr), cdrs) # count occurances per cdr3 for each motif

    end

    return motif_counts

end

# ✅ checked
# makes list of motifs in cdr3s, counts them, filters based on a cutoff, then provides set count
function get_motifs(st_::AbstractVector{<:AbstractString}, min::Int, max::Int)

    mo_ = Dict{String,Int}() # make dictionary

    @showprogress desc = "Finding and counting motifs..." for s1 in st_

        lastindex(s1) < 7 && continue # only keep cdr3 7 or longer

        s2 = s1[4:(end-3)] # remove first and last 3 aa
        for um in min:max

            um > lastindex(s2) && continue # make sure cdr3 is longer than the motif size

            i1 = 0
            i2 = i1 + um - 1

            while i2 < lastindex(s2)

                m = s2[(i1+=1):(i2+=1)] # get the motif
                mo_[m] = get!(mo_, m, 0) + 1 # count total motif occurances

            end

        end

    end

    mo_ = Dict(m => num for (m, num) in mo_ if num >= 3) # cutoff of 3
    mo2_ = get_motif_counts(collect(keys(mo_)), st_)
    return mo2_

end

function find_significant_motifs(motifs, cdrs1, cdrs2, nsim, ove_cutoff)

    Random.seed!(1)

    motifs_list = collect(keys(motifs))
    L = length(motifs_list)
    counts_orig = [motifs[m] for m in motifs_list]

    counts_sim = Array{Float64}(undef, nsim, length(motifs_list))

    significant_motifs = Dict{String,Float64}()

    @showprogress desc = "Calculating significant motifs..." for i in 1:nsim

        random_cdrs = sample(cdrs2, length(cdrs1); replace=true, ordered=false)
        random_counts = get_motif_counts(motifs_list, random_cdrs)
        for j in 1:L
            counts_sim[i, j] = get(random_counts, motifs_list[j], 0)
        end

    end

    # everything above this line works

    for (index, m) in enumerate(motifs_list)

        if ove_cutoff == true
            ove = counts_orig[index] / mean(counts_sim[:, index])

            if counts_orig[index] < 2
                continue

            elseif (counts_orig[index] == 2 && ove >= 1000) ||
                (counts_orig[index] == 3 && ove >= 100) ||
                (counts_orig[index] >= 4 && ove >= 10)

                wins = count(x -> x >= counts_orig[index], counts_sim[:, index])
                p_val = (wins + 1) / (nsim + 1)
                #print(p_val)
                if p_val <= 0.05
                    significant_motifs[m] = get!(significant_motifs, m, 0.0) + p_val
                end

            end
        end

        if ove_cutoff != true
            wins = count(x -> x >= counts_orig[index], counts_sim[:, index])
            p_val = (wins) / (nsim)
            print(p_val)
            if p_val <= 0.05
                significant_motifs[m] = get!(significant_motifs, m, 0.0) + p_val
            end
        end

    end

    println(string("Number of significant motifs identified:", length(significant_motifs)))
    return significant_motifs

end


#_________#

# takes list of cdrs and checks for a motif, making pairs
function make_motif_pairs(st_::AbstractVector{<:AbstractString}, motif)
    u1 = lastindex(st_)
    u2 = div(u1 * (u1 - 1), 2)

    in__ = Vector{Tuple{Int,Int}}(undef, u2)
    po_ = Vector{Int}(undef, u2)

    i1 = 0

    for i2 in 1:u1, i3 in (i2+1):u1
        s1, s2 = st_[i2], st_[i3]

        has1 = occursin(motif, s1)
        has2 = occursin(motif, s2)

        if has1 && has2
            i1 += 1
            in__[i1] = (i2, i3)
            po_[i1] = 1
        end
    end

    resize!(in__, i1)
    resize!(po_, i1)

    return in__, po_
end

#_________#

function get_motif(s1::AbstractString, um)

    s2 = s1[4:(end-3)]

    # TODO: Use Set
    st_ = String[]

    i1 = 0

    i2 = i1 + um - 1

    while i2 < lastindex(s2)

        push!(st_, s2[(i1+=1):(i2+=1)])

    end

    st_

end

function get_motif(st__, u1)

    di = Dict{String,Int}()

    for nd in eachindex(st__), st in st__[nd]

        if !haskey(di, st)

            di[st] = 0

        end

        di[st] += 1

    end

    [st for (st, u2) in di if u1 <= u2]

end


# =============================================================================================== #
# Graphing
# =============================================================================================== #

# TO DO: change label of vertices to the cell barcode, to allow for pairing analyses.

function make_cdr3_vertices!(g, cdrs)
    isblank(x) = x === nothing || x === missing || x in ("None", "none", "NA", "na", "Na", " ", "  ", "")

    for row in 1:nrow(cdrs)

        label = String(cdrs.cdr3[row])
        add_vertex!(g, label, Dict(:index => row)) # add vertex

        # check for vgenes, etc. and store if they exist
        if !isblank(cdrs.v_gene[row])
            g[label][:vgene] = cdrs.v_gene[row]
        end

        if !isblank(cdrs.j_gene[row])
            g[label][:jgene] = cdrs.j_gene[row]
        end

        if !isblank(cdrs.d_gene[row])
            g[label][:dgene] = cdrs.d_gene[row]
        end

    end

    return g # return graph

end

function make_barcode_vertices!(g, cdrs)
    isblank(x) = x === nothing || x === missing || x in ("None", "none", "NA", "na", "Na", " ", "  ", "", "   ")

    if :row_index ∉ names(cdrs)
            cdrs.row_index = collect(1:nrow(cdrs))
        end

    grouped_by_barcode = groupby(cdrs, :barcode)

    chain_rows = Dict{String, Vector{Int}}()

    for (i, subgroup) in enumerate(grouped_by_barcode)

        bc = String(sub.subgroup[1]) # get barcode for group

        if !has_vertex(g, bc)
            add_vertex!(g, bc, Dict(:index => i, :barcode => bc))
        end

        for row in 1:nrow(subgroup) # make dictionary for each vertex (barcode) of all chains for that barcode w/ row index

            if !isblank(subgroup.chain[row])
                chain = String(subgroup.chain[row])
            end

            chain_index = subgroup.row_index[row]  # original row index in `cdrs`

            chain_row_pair = get!(chain_rows, chain, Int[])
            push!(chain_row_pair, chain_index)

        end

        g[bc][:chains] = chain_rows # add dictionary to vertex

        if :sample ∈ names(subgroup) && !isblank(subgroup.sample[1]) # add sample if it exisrts
            g[bc][:sample] = subgroup.sample[1]
        end

    end # end of subgroup loop

    return g

end

# make a global edge w/ annotation of hamming disrance
function add_global_edge!(g, u, v, distance)
    add_edge!(g, u, v, Dict(:distance => distance))
end

# add a local edge w/ annotation of motif and the pval (from get_significant_motifs)
function add_local_edge!(g, u, v, motif::String, pval)
    add_edge!(g, u, v, Dict(:motif => motif, :mpval => pval))
end

# go through cdrs
function make_edges(cdrs, motifs, isglobal, islocal)

    cdr3_vec = String.(cdrs.cdr3)

    g = MetaGraph(
        Graph();
        label_type=String,
        vertex_data_type=Dict{Symbol,Any},  # test storing multiple labels?
        edge_data_type=Dict{Symbol,Any},
    )

    g = make_vertices!(g, cdrs)

    # start w/ global distances
    # # changed so label is cdr3
    if isglobal == true
        pairs, dists = make_distance(cdr3_vec)

            @showprogress desc = "Making edges..." for (index, (i, j)) in enumerate(pairs)
                d = dists[index]
                if d <= 1
                    # i, j are indices into cdr3_vec
                    u = cdr3_vec[i]
                    v = cdr3_vec[j]

                    if haskey(g, u, v)
                        g[u, v][:distance] = d
                    else
                        add_global_edge!(g, u, v, d)
                    end
                end
            end
    end

    # TODO: next do local edges
    if islocal == true

       @showprogress desc = "Making edges..." for (motif, pval) in motifs

            mask  = occursin.(Ref(motif), cdr3_vec)
            matches  = findall(mask)
            nmatches = length(matches)

            # nothing to do if fewer than 2 sequences have this motif
            if nmatches <= 1
                continue
            end

            # all unique pairs of those indices
            for i in 1:(nmatches-1), j in (i+1):nmatches
                u = cdr3_vec[matches[i]]
                v = cdr3_vec[matches[j]]

                if haskey(g, u, v)
                    g[u, v][:motif] = motif
                    g[u, v][:mpval] = pval
                else
                    add_local_edge!(g, u, v, motif, pval)
                end
            end
        end
    end # end of local edges

    return g

end

function get_cdr3s_for_vertex_chain(g, cdrs, bc::String, chain::String) # get cdr3s from graph indices for each barcode
    vertex = g[bc]
    haskey(vertex, :cdr3_rows_by_chain) || return String[]
    rows = get(vertex[:cdr3_rows_by_chain], chain, Int[])
    return String.(cdrs.cdr3[rows]) # retrieve list of all cdr3s for that barcode
end

function make_barcode_graph(g, cdrs)

    cdr3_vec = String.(cdrs.cdr3)

    g = MetaGraph(
        Graph();
        label_type=String,
        vertex_data_type=Dict{Symbol,Any},
        edge_data_type=Dict{Symbol,Any},
    )

    g = make_barcode_vertices!(g, cdrs)

    return g

end

function make_global_edges(g, cdrs; chain::String = "TRB", max_distance::Int = 1,)

    # subset dataframe to specific chain and get list of barcodes

    filtered_cdrs = subset(cdrs, :chain => ByRow(x -> x == chain))

    #barcodes = unique(filtered_cdrs.barcode)

    # using filtered df, make distances
    # go through index pairs and distances and check cutoff
    # if cutoff is good, find the row_index and add a global edge annotated with the specific chain

    cdr3s = String.(filtered_cdrs.cdr3)

    index_pairs, dists = make_distance(cdr3s)

    for (index, (i, j)) in enumerate(pairs)

        d = dists[index]

        if d <= 1

            bc_1 = String(filtered_cdrs.barcode[i])
            bc_2 = String(filtered_cdrs.barcode[j])
            cdr3_1 = String(filtered_cdrs.row_index[i])
            cdr3_2 = String(filtered_cdrs.row_index[j])

            bc_1 == bc_2 && continue

            if !has_vertex(g, bc_1)
                add_vertex!(g, bc_1, Dict(:barcode => bc_1))
            end

            if !has_vertex(g, bc_2)
                add_vertex!(g, bc_2, Dict(:barcode => bc_2))
            end

            add_global_edge!(g, bc_1, bc_2, d)
            g[bc_1, bc_2][:chains] = (cdr3_1, cdr3_2) # for now track the indices so that later people can go back to see pairs 4 the edges if they wany
        end

    end


end

function make_local_edges(g, cdrs; chain::String = "TRB", max_distance::Int = 1,)

    filtered_cdrs = subset(cdrs, :chain => ByRow(x -> x == chain))

    barcodes = unique(filtered_cdrs.barcode)



    # subset dataframe to specific chain and get list of barcodes

    # go through each vertex with those barcodes
    # filter chain dictionary to only the specific chain
    # then go through the indexes , i , i+1
    # checkthrough every motif to see if they have any shared and add a local edge for that
    # add a global edge for those with matches with a note on the cdr3 indexes involved in the edge



end

# =============================================================================================== #
# Clusters / Significance
#
# TO DO: multiple comparisons correction?
# =============================================================================================== #

# length score
# for each group, resample 1000 groups of group size n and measure cdr3 length distribution.

using StatsBase: mode

# 📋 TO DO: change to use barcodes and use cdr table as a lookup to get cdr3 length
function find_length_pvals(g, cdrs2, sim_depth)
    clusters = connected_components(g)

    modes = Vector{Int}()
    props = Vector{Float64}()

    # get most common length for each cluster
    for cluster in clusters

        labels = label_for.(Ref(g), cluster) # should be fixed now!
        cluster_lengths = [length(s) for s in labels]
        m = mode(cluster_lengths)
        push!(modes, m)
        push!(props, count(==(m), cluster_lengths) / length(cluster))

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

function score_lengths(g, cdrs2)

    clusters = connected_components(g)

    p_vals = find_length_pvals(g, cdrs2)

    for (index, cluster) in enumerate(clusters)

        labels = label_for.(Ref(g), cluster)

        for label in labels
            g[label][:length_pval] = p_vals[index] # adding cluster pval to each vertex
        end

    end

    return p_vals

end

#_______________#

# v gene score
# hyper geomtric or parameteric when over 200 sequences in one group

using Distributions

# 📋 TO DO: change to use barcodes and use cdr table to find vgenes for each pair? cell? individual chain?
# decide what to do re: that soon - might make sense to just score within all TRB, then TRA, TRG, etc. etc.
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
            vgene = g[label_for(g, vertex)][:vgene] # should be fixed now!
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

            # TO DO = verify if below is correct

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

        labels = label_for.(Ref(g), cluster)

        for label in labels
            g[label][:vgene_pval] = cluster_pvals[index] # adding cluster pval to each vertex
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

# 📋 TO DO: check to make sure its fully working using some mock data w/ fake donor names and the test cdr3 dataset
# also change to barcode system as above for other functions
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

            sample = g[label_for(g, vertex)][:sample]

            row_index = findfirst(==(sample), hla_df.donor)
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
