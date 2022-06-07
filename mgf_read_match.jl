using DataFrames
using BenchmarkTools
using StatsBase



function loadmgf(fname)
    FIELDS = ["TITLE=", "RTINSECONDS=", "PEPMASS=", "CHARGE=", "SCANS="]
    
    function format_precursor!(spectrum)
        if occursin( " ", spectrum["PEPMASS"] )
            spectrum["PEPMASS"] = map( x -> tryparse(Float64, x), split( spectrum["PEPMASS"] , " " ) )
        else
            spectrum["PEPMASS"] = [ tryparse(Float64, spectrum["PEPMASS"] ) ]
        end
        
        local polarity_multiplier = 1
        if spectrum["CHARGE"][end] == '-'
            polarity_multiplier = -1
        end
        
        if isnumeric( spectrum["CHARGE"][end] ) == false
            spectrum["CHARGE"] = polarity_multiplier * tryparse(Int32, spectrum["CHARGE"][1:end-1] )
        else
            spectrum["CHARGE"] = tryparse(Int32, spectrum["CHARGE"])
        end
        return true
    end
    
    function msdata_to_df!(spectrum)
        spectrum["ms_data"] = DataFrame( spectrum["ms_data"] )
        return true
    end
    
    local spectrum
    spectra_array = []
    open(fname) do fh
        state = false
        for line in eachline(fh)
            
            if length(line) > 0 && isnumeric( line[1] ) && state == true
                num_line = map( x -> tryparse(Float64, x), split(line, " ") )
                push!(
                    spectrum["ms_data"]["m/z"],
                    num_line[1]
                )
                push!(
                    spectrum["ms_data"]["Intensity"],
                    num_line[2]
                )
            elseif occursin("BEGIN IONS", line)
                spectrum = Dict{String,Any}( [ "ms_data" => Dict(
                            [
                                "m/z"=> Array{Float64}(undef, 0),
                                "Intensity"=> Array{Float64}(undef, 0)
                            ]
                        ) ] )
                state = true
                
            elseif occursin("END IONS", line)
                #println(spectrum)
                if length( spectrum["ms_data"] ) > 0
                    format_precursor!(spectrum)
                    msdata_to_df!(spectrum)
                    push!(spectra_array, spectrum)
                end
                state = false
            else
                for field_name in FIELDS
                    if occursin(field_name, line) && state == true
                        spectrum[ field_name[1:end-1] ] = split(line, field_name)[2]
                    end
                end
            end
        end
    end
    spectra_array
end

fname = "Yeast_1000spectra.mgf"

res = loadmgf(fname)
length(res)

res[501]

res[501]["ms_data"]

@benchmark loadmgf(fname) samples=200

AA_DELTAS = Dict(
    'G'=> 57.02147, 'A'=> 71.03712, 'S'=> 87.03203, 'P'=> 97.05277, 'V'=> 99.06842,
    'T'=> 101.04768, "Ccam"=> 160.03065, "Cmes"=> 148.996912, "I/L"=> 113.08407,
    'N'=> 114.04293, 'D'=> 115.02695, 'Q'=> 128.05858, 'K'=> 128.09497, 'E'=> 129.0426,
    'M'=> 131.04049, "Mox"=> 147.0354, 'H'=> 137.05891, 'F'=> 147.06842, 'R'=> 156.10112,
    'Y'=> 163.06333, 'W'=> 186.07932
)

#Flatten the values from the dictionary
single_res_Δ = collect( values(AA_DELTAS) )
println( typeof(single_res_Δ) )
#Add doubly-charged and triply-charged mass Deltas (simply divide by 2 and 3)
single_res_Δ = vcat( single_res_Δ / 3, single_res_Δ / 2, single_res_Δ )
println( length( single_res_Δ ) )
single_res_Δ[1:5]

function matches(spectra, masses_to_match, rel_tolerance)
    res_dict = Dict([
        ( "Spectrum_idx", Vector{Int64}() ),
        ( "Exp_idx", Vector{Int64}() ),
        ( "Library_idx", Vector{Int64}() ),
        ( "Rel_error", Vector{Float64}() )
    ])
    min_theo_val = minimum(masses_to_match) * (1 - rel_tolerance)
   for (idx, s) in enumerate(spectra)
        exp_Δs = filter(
            x -> x > min_theo_val,
            s["ms_data"][!, "m/z"] .- s["ms_data"][!, "m/z"]'
        )
        rel_deltas_matrix = map(
            abs,
            (masses_to_match .- exp_Δs')
        ) * 2 ./ (masses_to_match .+ exp_Δs')
        matching_inds = findall(x -> isless(x, rel_tolerance), rel_deltas_matrix)
        num_matches = size(matching_inds)[1]
        if num_matches > 0
            append!( res_dict["Spectrum_idx"], fill(idx, num_matches) )
            append!( res_dict["Library_idx"], map(x -> x[1], matching_inds) )
            append!( res_dict["Exp_idx"], map(x -> x[2], matching_inds) )
            append!( res_dict["Rel_error"], rel_deltas_matrix[matching_inds] )
        end
    end
    DataFrame(res_dict)
end

m = matches(res, single_res_Δ, 1e-5)
size(m)

@benchmark matches(res, single_res_Δ, 1e-5) samples=50



function matchesfuse(spectra, masses_to_match, rel_tolerance)
    res_dict = Dict([
        ( "Spectrum_idx", Vector{Int64}() ),
        ( "Exp_idx", Vector{Int64}() ),
        ( "Library_idx", Vector{Int64}() ),
        ( "Rel_error", Vector{Float64}() )
    ])
    min_theo_val = minimum(masses_to_match) * (1 - rel_tolerance)
   for (idx, s) in enumerate(spectra)
        exp_Δs = filter(
            x -> x > min_theo_val,
            s["ms_data"][!, "m/z"] .- s["ms_data"][!, "m/z"]'
        )
        #Fuse the vectorized operations using the dot macro
        rel_deltas_matrix = @. map(abs, (masses_to_match - exp_Δs') ) * 2 / (masses_to_match + exp_Δs')
        matching_inds = findall(x -> isless(x, rel_tolerance), rel_deltas_matrix)
        num_matches = size(matching_inds)[1]
        if num_matches > 0
            append!( res_dict["Spectrum_idx"], fill(idx, num_matches) )
            append!( res_dict["Library_idx"], map(x -> x[1], matching_inds) )
            append!( res_dict["Exp_idx"], map(x -> x[2], matching_inds) )
            append!( res_dict["Rel_error"], rel_deltas_matrix[matching_inds] )
        end
    end
    DataFrame(res_dict)
end

m = matchesfuse(res, single_res_Δ, 1e-5)
size(m)

m[ m.Spectrum_idx .== 2 , :]

@benchmark matchesfuse(res, single_res_Δ, 1e-5) samples=50

aa_curated = [
    'G', 'A', 'S', 'P', 'V', 'T', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'
]

sample(aa_curated)

a = 'H'
a *= 'A'
a

a = ""
for i = 1:20
    a *= sample(aa_curated)
end
a

length(a)

sequences_list = Array([])
for i = 1:10000
    a = ""
    for j = 1:20
        a *= sample(aa_curated)
    end
    append!(sequences_list, [a])
end
length(sequences_list)

sequences_list[1:5]

for i in split("ERPPPGVMGDRSYPVRFVDP", "")
    println(AA_DELTAS[i])
end

function calculate_masses_loop(sequences_list)
    masses = Array{Float64}([])
    for i in sequences_list
        m = 18.010565
        for j in i
            m += AA_DELTAS[j]
        end
        append!(masses, [m])
    end
    masses
end

l = calculate_masses_loop(sequences_list)
length(l)

l[1:5]

@benchmark calculate_masses_loop(sequences_list) samples=200
