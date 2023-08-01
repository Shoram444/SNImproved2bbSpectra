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
    maxSampleSize = 1_000_000,
    verbose = true
) where T<:Real
    pVals = Vector{T}(undef, nSamples)   # initiaite a container to hold vectors of 100 p-values 
    efficiency = float(0)                # efficiency 
    nInARowCounter = 0                   # a counter for n in a row 100% efficiencies

    samples1 = Vector{T}(undef, nSamples)
    samples2 = Vector{T}(undef, nSamples)

    i = 1 
    sampleSize = sampleSizes[i]

    while (nInARowCounter != nInARow && sampleSize <= maxSampleSize)
        if( i <= length(sampleSizes) )
            sampleSize = sampleSizes[i]
        else
            sampleSize += 50_000 
        end

        samples1 = get_samples(vector1, sampleSize, nSamples, true)
        samples2 = get_samples(vector2, sampleSize, nSamples, true)

        @views pVals = ApproximateTwoSampleKSTest.(samples1, samples2) .|> pvalue
        efficiency = get_efficiency(pVals, CL, nSamples)

        if(efficiency == 1.0)
            nInARowCounter += 1
        else
            nInARowCounter = 0
        end

        if(verbose)
            @show "Current Sample Size = $sampleSize; ε = $efficiency"
        end

        i += 1
    end
    minEvents = sampleSize < maxSampleSize ? sampleSize : maxSampleSize

    if(verbose)
        if( sampleSize > maxSampleSize )
            @warn "Did not reach $nInARow 100% efficiencies before reaching max sample size= $maxSampleSize. Setting best sample size to $maxSampleSize."
        end
        @show "best sample size: $minEvents"
    end

    return KS(
        vector1,
        vector2,
        CL,
        minEvents
    )
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


function get_pVals_Fast(ks::KS, sampleSizes, nSamples = 100)
    pVals = Vector{Vector{<:Real}}(undef, length(sampleSizes))   # initiaite a container to hold vectors of 100 p-values for each sample size
    
    @inbounds Threads.@threads for i in eachindex(sampleSizes)

        pVals[i] = ApproximateTwoSampleKSTest.(
            get_samples(ks.vector1, sampleSizes[i], nSamples, true), 
            get_samples(ks.vector2, sampleSizes[i], nSamples, true)
            ) .|> pvalue
    end

    return pVals
end