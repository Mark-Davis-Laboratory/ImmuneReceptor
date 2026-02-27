# preparing gamma delta naive dataset

using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("StatsBase")
Pkg.add("Plots")
using CSV: read
using Plots
using DataFrames
using StatsBase

# TO-DO:
# remove cells without TRG, TRD or CDR3 read ✅
# grab TRJ, TRV, TRD and CDR3
# calculate percent utilization for each V, D, and J gene for G and D chains ✅
# calculate distribtuion of lengths for CDR3s ✅
# maybe make symbol/AA utilization plot for CDR3s?

# load datasets (load just one for now to test)
df = CSV.read(
    "/Users/dr/Downloads/GSM6470059_PNT3_filtered_contig_annotations_copy.csv",
    DataFrame,
)

# remove any rows that don't have TRG or TRD value in "chain" column or without a value in the CDR3 column

filtered_df = df[((df.chain.=="TRG").|(df.chain.=="TRD")).&(df.cdr3.!="None"), :]

# count incidence of every unique value in TRJ, TRV or TRD column.

df_G = filter(:chain => ==("TRG"), filtered_df)
df_G = filter(:v_gene => !=("None"), df_G)
df_G = filter(:j_gene => !=("None"), df_G)

df_D = filter(:chain => ==("TRD"), filtered_df)
df_D = filter(:v_gene => !=("None"), df_D)
df_D = filter(:j_gene => !=("None"), df_D)
df_D = filter(:d_gene => !=("None"), df_D)

vgenes_G = proportionmap(df_G.v_gene)
vgenes_D = proportionmap(df_D.v_gene)
jgenes_G = proportionmap(df_G.j_gene)
jgenes_D = proportionmap(df_D.j_gene)
dgenes_G = proportionmap(df_G.d_gene)
dgenes_D = proportionmap(df_D.d_gene)

# make plot

# γ-chain plots
p1 = bar(
    collect(keys(vgenes_G)),
    collect(values(vgenes_G));
    legend=false,
    xlabel="V Gene",
    ylabel="Frequency",
    title="γ-Chain V Gene Frequency",
)

p2 = bar(
    collect(keys(jgenes_G)),
    collect(values(jgenes_G));
    legend=false,
    xlabel="J Gene",
    ylabel="Frequency",
    title="γ-Chain J Gene Frequency",
)

# δ-chain plots
p3 = bar(
    collect(keys(vgenes_D)),
    collect(values(vgenes_D));
    legend=false,
    xlabel="V Gene",
    ylabel="Frequency",
    title="δ-Chain V Gene Frequency",
)

p4 = bar(
    collect(keys(dgenes_D)),
    collect(values(dgenes_D));
    legend=false,
    xlabel="D Gene",
    ylabel="Frequency",
    title="δ-Chain D Gene Frequency",
)

p5 = bar(
    collect(keys(jgenes_D)),
    collect(values(jgenes_D));
    legend=false,
    xlabel="J Gene",
    ylabel="Frequency",
    title="δ-Chain J Gene Frequency",
)

plot(p1, p2, p3, p4, p5; layout=(2, 3))

# count incidence of every unique length for unique CDR3s

lengths_D = countmap(length.(df_D.cdr3))
lengths_G = countmap(length.(df_G.cdr3))

p1_b = bar(
    collect(keys(lengths_D)),
    collect(values(lengths_D));
    legend=false,
    xlabel="CDR3 Length (in AA)",
    ylabel="Frequency",
    title="δ-Chain CDR3 Lengths",
)

p2_b = bar(
    collect(keys(lengths_G)),
    collect(values(lengths_G));
    legend=false,
    xlabel="CDR3 Length (in AA)",
    ylabel="Frequncy",
    title="γ-Chain CDR3 Lengths",
)

plot(p1_b, p2_b)

# make sequence logo for CDR3s (maybe center sequences)
# TO DO: figure out how to make sequene logo with centered sequences

function center_cdr3s(s::AbstractString, maximum_length::Int; padchar::Char='-')
    len = length(s)
    padding_left = div(maximum_length - len, 2)
    padding_right = maximum_length - len - padding_left
    return repeat(string(padchar), padding_left) *
           s *
           repeat(string(padchar), padding_right)
end

function center_cdr3s(strings::Vector{<:AbstractString}; padchar::Char='-')
    valid_cdrs = filter(x -> x != "none", strings)
    if isempty(valid_cdrs)
        print("No valid CDRs")
        return String[]
    end

    maximum_length = maximum(length.(valid_cdrs))

    vals = Vector{String}(undef, length(valid_cdrs))
    for (i, cdr3) in enumerate(valid_cdrs)
        vals[i] = center_cdr3s(cdr3, maximum_length; padchar=padchar)
    end
    return vals
end

centered_d = center_cdr3s(df_D.cdr3)
centered_g = center_cdr3s(df_G.cdr3)
