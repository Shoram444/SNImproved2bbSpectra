function get_min_value(mat::Matrix{<:Real})
    flatArray = collect(Iterators.flatten(mat))

    return return (minimum(flatArray[.!isnan.(flatArray)]))::Real
end

function get_min_idx(mat::Matrix{<:Real})
    minVal = Inf
    idx = (0, 0)
    for i in CartesianIndices(mat)
        if (mat[i] != NaN && mat[i] < minVal)
            minVal = mat[i]
            idx = i
        end
    end
    return idx::CartesianIndex{2}
end

function cartesianIdx_to_range(idx::CartesianIndex{2}, Δϕ::Real)
    range::Tuple{<:Real,<:Real} = idx[2] * Δϕ - Δϕ, idx[1] * Δϕ
    return range
end

function get_best_sample_size(
    efficiencies::Vector{<:Real},
    sampleSizes::Vector{<:Real},
    nInARow::Int = 3,
)
    idx = 1

    length(efficiencies) < nInARow &&
        error("size of data: $(length(efficiencies)) is less than nInARow: $nInARow")

    sumOfLastN = sum(efficiencies[end-nInARow+1:end]) # sum of the last nInARow numbers

    if (sumOfLastN < nInARow)
        @warn(
            "Did not get $nInARow 100% efficiencies in a row. Please increase sample size. Returning last sampleSize."
        )
        return sampleSizes[end]
    end

    for j = length(efficiencies):-1:nInARow+1 # iterate backwards
        sumOfEffs = sum(efficiencies[j-nInARow+1:j])
        idx = j + 1
        if (sumOfEffs < nInARow)
            break
        end
    end
    return sampleSizes[idx]
end


@recipe function f(h::Histogram{T,1,Tuple{Vector{Int64}}}) where {T<:Real}
    bins = vcat([h.edges[1][1]], midpoints(h.edges[1]), [h.edges[1][end]])

    bincounts = vcat([h.weights[1]], h.weights, [h.weights[end]])

    seriestype --> :stepmid
    x := bins
    y := bincounts
    ()
end


function fill_from_root_file(inFile::ROOTFile, treeName::String, key::String)
    data = Float64[]

    if (key ∉ keys(inFile[treeName]))
        error(" Key: $key is not in tree: $tree. Provide a valid key! ")
    end

    lt = LazyTree(inFile, treeName, [key])

    for evt in lt
        push!(data, evt[1])
    end
    return data
end

function get_samples(
    data::Vector{<:Real},
    sample_size::Real,
    n_samples::Real,
    replace = false,
)
    sample_size = Int(sample_size)
    n_samples = Int(n_samples)
    data_size = length(data)

    if (replace == false && sample_size * n_samples > data_size)
        error(
            "Too many samples, $sample_size*$n_samples is larger than size of data: $data_size",
        )
    end

    subsets = Array{Vector{Int}}(undef, n_samples)
    if (replace == false)
        sample_indices = StatsBase.sample(1:data_size, sample_size, replace = replace)
        subsets = [data[sample_size*(x-1)+1:sample_size*x] for x = 1:n_samples]
    else
        subsets =
            [StatsBase.sample(data, sample_size, replace = replace) for x = 1:n_samples]
    end

    return subsets
end

function get_efficiency(vector::Vector{<:Real}, CL::Real, nSamples = length(vector))
    return count(p -> p <= 1 - CL, vector) / nSamples
end

function get_residuals(vector1, vector2, binrange, normed = true)
    if (normed)
        bincounts1 = StatsBase.fit(Histogram{Float64}, vector1, binrange) |> normalize
        bincounts2 = StatsBase.fit(Histogram{Float64}, vector2, binrange) |> normalize
    else
        bincounts1 = StatsBase.fit(Histogram{Float64}, vector1, binrange)
        bincounts2 = StatsBase.fit(Histogram{Float64}, vector2, binrange)
    end
    bincounts1.weights
    bincounts2.weights

    return bincounts2.weights ./ bincounts1.weights
end


function get_residuals_errors(vector1, vector2, binrange, normed = false )
	hist1 = StatsBase.fit(Histogram{Float64}, vector1, binrange) 
	hist2 = StatsBase.fit(Histogram{Float64}, vector2, binrange) 

	bincounts1 = hist1.weights
	bincounts2 = hist2.weights
	
	if (normed)
		nTotal1 = length(vector1)
        nTotal2 = length(vector2)

		bincounts1 ./= nTotal1 # normalized bincounts
        bincounts2 ./= nTotal2 # normalized bincounts
	end
	r = bincounts1./ bincounts2
	res = @. r*sqrt( 1/ bincounts1 + 1/ bincounts2)

	return res
end

function prop_error(npassed, ntot, percent = false)
    if (ntot != 100)
        error("Uncertainties were calculated only for N=100, provided N=$ntot!")
    end

    errors100 = [
        0.000000e+00 0.000000e+00 1.121807e-02
        1.000000e-02 2.785118e-03 2.456629e-02
        2.000000e-02 8.872022e-03 3.791181e-02
        3.000000e-02 1.592492e-02 5.066677e-02
        4.000000e-02 2.334868e-02 6.289157e-02
        5.000000e-02 3.112133e-02 7.485862e-02
        6.000000e-02 3.921413e-02 8.669834e-02
        7.000000e-02 4.744811e-02 9.832986e-02
        8.000000e-02 5.580855e-02 1.098041e-01
        9.000000e-02 6.430260e-02 1.211748e-01
        1.000000e-01 7.280360e-02 1.323494e-01
        1.100000e-01 8.157604e-02 1.436188e-01
        1.200000e-01 9.018141e-02 1.545652e-01
        1.300000e-01 9.897695e-02 1.655619e-01
        1.400000e-01 1.079330e-01 1.765921e-01
        1.500000e-01 1.169662e-01 1.875845e-01
        1.600000e-01 1.259758e-01 1.984475e-01
        1.700000e-01 1.349889e-01 2.092161e-01
        1.800000e-01 1.441486e-01 2.200400e-01
        1.900000e-01 1.532961e-01 2.307665e-01
        2.000000e-01 1.624585e-01 2.414277e-01
        2.100000e-01 1.718053e-01 2.521980e-01
        2.200000e-01 1.808864e-01 2.626301e-01
        2.300000e-01 1.900989e-01 2.731261e-01
        2.400000e-01 1.995199e-01 2.837646e-01
        2.500000e-01 2.089822e-01 2.943826e-01
        2.600000e-01 2.181381e-01 3.046348e-01
        2.700000e-01 2.277205e-01 3.152546e-01
        2.800000e-01 2.372144e-01 3.257311e-01
        2.900000e-01 2.466333e-01 3.360787e-01
        3.000000e-01 2.561578e-01 3.464801e-01
        3.100000e-01 2.656193e-01 3.567676e-01
        3.200000e-01 2.752097e-01 3.671350e-01
        3.300000e-01 2.847895e-01 3.774437e-01
        3.400000e-01 2.943707e-01 3.877069e-01
        3.500000e-01 3.039832e-01 3.979556e-01
        3.600000e-01 3.136907e-01 4.082545e-01
        3.700000e-01 3.232659e-01 4.183768e-01
        3.800000e-01 3.330313e-01 4.286461e-01
        3.900000e-01 3.427805e-01 4.388567e-01
        4.000000e-01 3.523695e-01 4.488647e-01
        4.100000e-01 3.620832e-01 4.589562e-01
        4.200000e-01 3.819674e-01 4.817691e-01
        4.300000e-01 3.817180e-01 4.792238e-01
        4.400000e-01 3.915312e-01 4.892929e-01
        4.500000e-01 4.012813e-01 4.992590e-01
        4.600000e-01 4.111368e-01 5.092910e-01
        4.700000e-01 4.210735e-01 5.193646e-01
        4.800000e-01 4.308871e-01 5.292759e-01
        4.900000e-01 4.407896e-01 5.392371e-01
        5.000000e-01 4.507300e-01 5.491970e-01
        5.100000e-01 4.607081e-01 5.591555e-01
        5.200000e-01 4.706580e-01 5.690469e-01
        5.300000e-01 4.806541e-01 5.789452e-01
        5.400000e-01 4.906886e-01 5.888428e-01
        5.500000e-01 5.006331e-01 5.986110e-01
        5.600000e-01 5.107426e-01 6.085044e-01
        5.700000e-01 5.207335e-01 6.182394e-01
        5.800000e-01 5.308928e-01 6.281026e-01
        5.900000e-01 5.410829e-01 6.379562e-01
        6.000000e-01 5.509821e-01 6.474777e-01
        6.100000e-01 5.729496e-01 6.728868e-01
        6.200000e-01 5.714176e-01 6.670326e-01
        6.300000e-01 5.816934e-01 6.768047e-01
        6.400000e-01 5.918487e-01 6.864127e-01
        6.500000e-01 6.020150e-01 6.959877e-01
        6.600000e-01 6.121861e-01 7.055228e-01
        6.700000e-01 6.225794e-01 7.152338e-01
        6.800000e-01 6.328929e-01 7.248185e-01
        6.900000e-01 6.431387e-01 7.342874e-01
        7.000000e-01 6.535195e-01 7.438421e-01
        7.100000e-01 6.639322e-01 7.533781e-01
        7.200000e-01 6.742636e-01 7.627807e-01
        7.300000e-01 6.847801e-01 7.723147e-01
        7.400000e-01 6.952072e-01 7.817037e-01
        7.500000e-01 7.056463e-01 7.910470e-01
        7.600000e-01 7.162040e-01 8.004492e-01
        7.700000e-01 7.266659e-01 8.096932e-01
        7.800000e-01 7.372987e-01 8.190429e-01
        7.900000e-01 7.478892e-01 8.282821e-01
        8.000000e-01 7.585576e-01 8.375274e-01
        8.100000e-01 7.692570e-01 8.467280e-01
        8.200000e-01 7.799368e-01 8.558289e-01
        8.300000e-01 7.906863e-01 8.649141e-01
        8.400000e-01 8.015840e-01 8.740563e-01
        8.500000e-01 8.124342e-01 8.830531e-01
        8.600000e-01 8.233340e-01 8.919935e-01
        8.700000e-01 8.342813e-01 9.008663e-01
        8.800000e-01 8.453380e-01 9.097222e-01
        8.900000e-01 8.564406e-01 9.184840e-01
        9.000000e-01 8.675766e-01 9.271231e-01
        9.100000e-01 8.788755e-01 9.357484e-01
        9.200000e-01 8.902361e-01 9.442324e-01
        9.300000e-01 9.017385e-01 9.526210e-01
        9.400000e-01 9.133407e-01 9.608256e-01
        9.500000e-01 9.251661e-01 9.689042e-01
        9.600000e-01 9.371470e-01 9.766907e-01
        9.700000e-01 9.494374e-01 9.841794e-01
        9.800000e-01 9.621476e-01 9.911880e-01
        9.900000e-01 9.754569e-01 9.972389e-01
        1.000000e+00 9.887819e-01 1.000000e+00
    ]

    lower = npassed / 100 - errors100[Int(npassed)+1, 2]
    upper = errors100[Int(npassed)+1, 3] - npassed / 100

    return percent ? (lower * 100, upper * 100) : (lower, upper)
end


const G0 = 0.334863e-46 # [MeV]
const G2 = 0.148350E-46 # [MeV]
const G22 = 0.188606E-47 # [MeV]
const G4 = 0.824467E-47 # [MeV]

const H0 = 0.226022E-46 # [MeV]
const H2 = 0.929409E-47 # [MeV]
const H22 = 0.108907E-47 # [MeV]
const H4 = 0.484671E-47 # [MeV]

const ξ51 = 0.1397 # for SSD


function get_kappa(
    ξ31,
    ξ51 = 0.1397,
    G0 = 0.334863e-46,
    G2 = 0.148350E-46,
    G22 = 0.188606E-47,
    G4 = 0.824467E-47,
    H0 = 0.226022E-46,
    H2 = 0.929409E-47,
    H22 = 0.108907E-47,
    H4 = 0.484671E-47,
)
    numerator = H0 + ξ31 * H2 + 5 / 9 * ξ31^2 * H22 + (2 / 9 * ξ31^2 + ξ51) * H4
    denumerator = G0 + ξ31 * G2 + 1 / 3 * ξ31^2 * G22 + (1 / 3 * ξ31^2 + ξ51) * G4

    return -1 * numerator / denumerator
end

function DeltaM(r, M, N)
    sq = sqrt(1 / M + 1 / N)
    num = 3r + sqrt(r) * (r + 3) + 1
    den = (1 - r)^3

    return num / den * r * sq
end

function Mmin(r, M, N)
    num = r + sqrt(r) * sqrt(M / N)
    den = 1 - r

    return (num / den)^2
end
