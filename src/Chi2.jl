import Distributions: Chisq, ccdf


mutable struct Chi2{T}
	vector1::Vector{T}
	vector2::Vector{T}
	xMin::Real # minimum value of the array (i.e. for angles it would be 0 degrees)
	xMax::Real # maximum value of the array (i.e. for angles it would be 180 degrees)
	stepSize::Real
	sampleSizes::Vector{<:Real}
	CL::Real
	pVals::Vector{Vector{<:Real}}
	efficiencies::Vector{<:Real}
	minEvents::Real

	function Chi2(vector1::Vector{T}, vector2::Vector{T}, xMin::Real, xMax::Real, stepSize::Real, sampleSizes::Vector{<:Real}, CL::Real) where T
		if( (xMax-xMin)/stepSize%1 != 0.0 )  						# check the boundaries of the histogram
			error("Range must be integer divisible by stepsize! ")
		end
		nSamples = 100

		pVals = Vector{Vector{<:Real}}(undef, length(sampleSizes))   # initiaite a container to hold vectors of 100 p-values for each sample size
		efficiencies = Vector{Real}(undef, length(sampleSizes))	 # initiaite a container to hold the efficiency for each sample size

		for (i,sampleSize) in enumerate(sampleSizes)
			samples1 = get_samples(vector1, sampleSize, nSamples, true)
			samples2 = get_samples(vector2, sampleSize, nSamples, true)

			pVals[i] = ChisqTest.(samples1, samples2, xMin, xMax, stepSize) .|> pvalue
			efficiencies[i] = get_efficiency(pVals[i], CL, nSamples) 
		end

		minEvents = get_best_sample_size(efficiencies, sampleSizes, 3)

		new{T}(
				vector1::Vector{T}, 
				vector2::Vector{T}, 
				xMin::Real, 
				xMax::Real, 
				stepSize::Real, 
				sampleSizes::Vector{<:Real},
				CL::Real, 	
				pVals::Vector{Vector{<:Real}}, 
				efficiencies::Vector{<:Real}, 
				minEvents::Real
			)
	end
end

"""
ChisqTest(vector1, vector2, xMin, xMax, nBins)

A modified version of HypothesisTests.ChisqTest() which takes as an input 2 vectors of measured data: vector1, vector2 and compares their distributions in a histogram form.

Input parameters:
* vector1, vector2 are vectors of measured values
* xMin, xMax, binStep specify the minimum, maximum and step of the bins to compare the distributions

"""
function ChisqTest(vector1::Vector{<:Real}, vector2::Vector{<:Real}, xMin::Real, xMax::Real, stepSize::Real) 
    if( (xMax-xMin)/stepSize%1 != 0.0 )  						# check the boundaries of the histogram
        error("Range must be integer divisible by stepsize! ")
    end
	
	binRange = xMin:stepSize:xMax
	
	bincounts1 = bincounts(Hist1D(vector1, (binRange))) 
	bincounts2 = bincounts(Hist1D(vector2, (binRange))) 

    (minimum(bincounts1) < 5 || minimum(bincounts1) <5 )  && error("ERROR: method only works for bins with >5 counts, please edit histogram range.")
	
	dof = length(bincounts1) - 1  # degrees of freedom

    if ( length(vector1) == length(vector2))
        testStatistic = sum((bincounts2 - bincounts1).^2 ./bincounts1)
    else

        n1 = sum(bincounts1)
        n2 = sum(bincounts2)
        
        p1 = bincounts1 ./ n1 # get the normalized bincounts expectation
        p2 = bincounts2 ./ n2 
	
        testStatistic = sum((p2 - p1).^2 ./p1)
    end

    pvalue = ccdf( Chisq(dof), testStatistic )
	return testStatistic, pvalue, dof
end


function HypothesisTests.pvalue(ChiTest::Tuple{Real, Real, Real})
	return ChiTest[2]
end


function get_best_sample_size(efficiencies::Vector{<:Real}, sampleSizes::Vector{<:Real}, nInARow::Int = 3)
	idx = 1
	
	length(efficiencies) < nInARow && error("size of data: $(length(efficiencies)) is less than nInARow: $nInARow")

	sumOfLastN = sum(efficiencies[ end-nInARow+1:end]) # sum of the last nInARow numbers

	if (sumOfLastN < nInARow)
		@warn("Did not get $nInARow 100% efficiencies in a row. Please increase sample size. Returning last sampleSize.")
		return sampleSizes[end]
	end

	for j in length(efficiencies):-1:nInARow+1 # iterate backwards
		sumOfEffs = sum(efficiencies[j-nInARow+1:j])
        idx = j+1
		if ( sumOfEffs < nInARow )
            break
		end
	end
	return sampleSizes[idx]
end

function AnalysisModule.get_efficiency(chi2::Chi2, CL::Real, sampleSize::Real) 
	pValsId = isempty(findfirst( x -> x == sampleSize, chi2.sampleSizes )) ? error("invalid sample size, $sampleSize, not in provided sample sizes $chi2.sampleSizes") : findfirst( x -> x == sampleSize, chi2.sampleSizes ) 

	return AnalysisModule.get_efficiency(chi2.pVals[pValsId], chi2.CL) 
end

function get_best_sample_size(chi::Chi2, sampleSizes, CL, nSamples, nInARow)
	effs = get_efficiency.(chi.pVals, CL, nSamples) 
	return get_best_sample_size(effs, sampleSizes, nInARow)
end
