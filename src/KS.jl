mutable struct KS
    vector1::Vector{<:Real}
    vector2::Vector{<:Real}
    CL::Real
    minEvents::Real
end

function KS(
    vector1::Vector{T},
    vector2::Vector{T},
    CL::Real;
    nInARow = 3,
    nSamples = 100,
    sampleSizes = vcat(                 #default sample sizes, if different values are required, they can be passed as a new vector
        collect(10_000:1000:19_000),
        collect(20_000:10_000:100_000),
        collect(100_000:50_000:950_000),
    ),
    maxSampleSize = 1_000_000
) where T<:Real
    pVals = Vector{T}(undef, nSamples)   # initiaite a container to hold vectors of 100 p-values 
    efficiency = float(0)                # efficiency 
    nInARowCounter = 0                   # a counter for n in a row 100% efficiencies

    samples1 = Vector{T}(undef, nSamples)
    samples2 = Vector{T}(undef, nSamples)

    i = 1 
    sampleSize = sampleSizes[i]

    while (nInARowCounter != nInARow || sampleSize == maxSampleSize)
        sampleSize = sampleSizes[i]
        samples1 = get_samples(vector1, sampleSize, nSamples, true)
        samples2 = get_samples(vector2, sampleSize, nSamples, true)

        @views pVals = ApproximateTwoSampleKSTest.(samples1, samples2) .|> pvalue
        efficiency = get_efficiency(pVals, CL, nSamples)

        if(efficiency == 1.0)
            nInARowCounter += 1
        else
            nInARowCounter = 0
        end
        i += 1
    end
    minEvents = sampleSize

    return KS(
        vector1,
        vector2,
        CL,
        minEvents
    )
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


function get_pVals(ks::KS, sampleSizes, nSamples = 100)
    pVals = Vector{Vector{<:Real}}(undef, length(sampleSizes))   # initiaite a container to hold vectors of 100 p-values for each sample size
    
    for (i, sampleSize) in enumerate(sampleSizes)
        samples1 = get_samples(ks.vector1, sampleSize, nSamples, true)
        samples2 = get_samples(ks.vector2, sampleSize, nSamples, true)

        pVals[i] = ApproximateTwoSampleKSTest.(samples1, samples2) .|> pvalue
    end

    return pVals
end

# function AnalysisModule.get_efficiency(ks::KS, CL::Real, sampleSize::Real)
#     pValsId =
#         isempty(findfirst(x -> x == sampleSize, ks.sampleSizes)) ?
#         error(
#             "invalid sample size, $sampleSize, not in provided sample sizes $ks.sampleSizes",
#         ) : findfirst(x -> x == sampleSize, chi2.sampleSizes)

#     return AnalysisModule.get_efficiency(ks.pVals[pValsId], ks.CL)
# end

# function get_best_sample_size(ks::KS, sampleSizes, CL, nSamples, nInARow)
#     effs = get_efficiency.(ks.pVals, CL, nSamples)
#     return get_best_sample_size(effs, sampleSizes, nInARow)
# end
