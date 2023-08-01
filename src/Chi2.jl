import Distributions: Chisq, ccdf


mutable struct Chi2{T}
    vector1::Vector{T}
    vector2::Vector{T}
    CL::Real
    minEvents::Real
    xMin::Real
    xMax::Real
    xStep::Real
end

function Chi2(
    vector1::Vector{T},
    vector2::Vector{T},
    CL::Real,
    xMin::Real,
    xMax::Real,
    xStep::Real;
    nInARow = 3,
    nSamples = 100,
    sampleSizes = vcat(                 #default sample sizes, if different values are required, they can be passed as a new vector
        collect(10_000:1000:19_000),
        collect(20_000:10_000:100_000),
        collect(100_000:50_000:950_000),
    ),
    maxSampleSize = 1_000_000,
    verbose = :true
) where T<:Real
    if ((xMax - xMin) / xStep % 1 != 0.0)  # check the boundaries of the histogram
        error("Range must be integer divisible by xStep! ")
    end

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

        @views pVals = ChisqTest.(samples1, samples2, xMin, xMax, xStep) .|> pvalue
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
    
    return Chi2(
        vector1,
        vector2,
        CL,
        minEvents,
        xMin,
        xMax,
        xStep
    )
end



"""
ChisqTest(vector1, vector2, xMin, xMax, nBins)

A modified version of HypothesisTests.ChisqTest() which takes as an input 2 vectors of measured data: vector1, vector2 and compares their distributions in a histogram form.

Input parameters:
* vector1, vector2 are vectors of measured values
* xMin, xMax, binStep specify the minimum, maximum and step of the bins to compare the distributions

"""
function ChisqTest(
    vector1::Vector{<:Real},
    vector2::Vector{<:Real},
    xMin::Real,
    xMax::Real,
    xStep::Real,
)
    if ((xMax - xMin) / xStep % 1 != 0.0)  # check the boundaries of the histogram
        error("Range must be integer divisible by xStep! ")
    end

    binRange = xMin:xStep:xMax

    bincounts1 = bincounts(Hist1D(vector1, (binRange)))
    bincounts2 = bincounts(Hist1D(vector2, (binRange)))

    (minimum(bincounts1) < 5 || minimum(bincounts1) < 5) && error(
        "ERROR: method only works for bins with >5 counts, please edit histogram range.",
    )

    dof = length(bincounts1) - 1  # degrees of freedom

    if (length(vector1) == length(vector2))
        testStatistic = sum((bincounts2 - bincounts1) .^ 2 ./ bincounts1)
    else

        n1 = sum(bincounts1)
        n2 = sum(bincounts2)

        p1 = bincounts1 ./ n1 # get the normalized bincounts expectation
        p2 = bincounts2 ./ n2

        testStatistic = sum((p2 - p1) .^ 2 ./ p1)
    end

    pvalue = ccdf(Chisq(dof), testStatistic)
    return testStatistic, pvalue, dof
end


function HypothesisTests.pvalue(ChiTest::Tuple{Real,Real,Real})
    return ChiTest[2]
end



function AnalysisModule.get_efficiency(chi2::Chi2, CL::Real, sampleSize::Real)
    pValsId =
        isempty(findfirst(x -> x == sampleSize, chi2.sampleSizes)) ?
        error(
            "invalid sample size, $sampleSize, not in provided sample sizes $chi2.sampleSizes",
        ) : findfirst(x -> x == sampleSize, chi2.sampleSizes)

    return AnalysisModule.get_efficiency(chi2.pVals[pValsId], chi2.CL)
end

function get_best_sample_size(chi::Chi2, sampleSizes, CL, nSamples, nInARow)
    effs = get_efficiency.(chi.pVals, CL, nSamples)
    return get_best_sample_size(effs, sampleSizes, nInARow)
end


function get_pVals(chi2::Chi2, sampleSizes, nSamples = 100)
    pVals = Vector{Vector{<:Real}}(undef, length(sampleSizes))   # initiaite a container to hold vectors of 100 p-values for each sample size
    
    for (i, sampleSize) in enumerate(sampleSizes)
        samples1 = get_samples(chi2.vector1, sampleSize, nSamples, true)
        samples2 = get_samples(chi2.vector2, sampleSize, nSamples, true)

        pVals[i] = ChisqTest.(samples1, samples2, chi2.xMin, chi2.xMax, chi2.xStep) .|> pvalue
    end

    return pVals
end

function get_pVals_Fast(chi2::Chi2, sampleSizes, nSamples = 100)
    pVals = Vector{Vector{<:Real}}(undef, length(sampleSizes))   # initiaite a container to hold vectors of 100 p-values for each sample size
    
    @inbounds Threads.@threads for i in eachindex(sampleSizes)

        pVals[i] = ChisqTest.(
            get_samples(chi2.vector1, sampleSizes[i], nSamples, true), 
            get_samples(chi2.vector2, sampleSizes[i], nSamples, true), 
            chi2.xMin, 
            chi2.xMax, 
            chi2.xStep
            ) .|> pvalue
    end

    return pVals
end