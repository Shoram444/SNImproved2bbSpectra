mutable struct KS{T}
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

	function KS(vector1::Vector{T}, vector2::Vector{T}, xMin::Real, xMax::Real, stepSize::Real, sampleSizes::Vector{<:Real}, CL::Real) where T
		if( (xMax-xMin)/stepSize%1 != 0.0 )  						# check the boundaries of the histogram
			error("Range must be integer divisible by stepsize! ")
		end
		nSamples = 100

		pVals = Vector{Vector{<:Real}}(undef, length(sampleSizes))   # initiaite a container to hold vectors of 100 p-values for each sample size
		efficiencies = Vector{Real}(undef, length(sampleSizes))	 # initiaite a container to hold the efficiency for each sample size

		for (i,sampleSize) in enumerate(sampleSizes)
			samples1 = get_samples(vector1, sampleSize, nSamples, true)
			samples2 = get_samples(vector2, sampleSize, nSamples, true)

			pVals[i] = ApproximateTwoSampleKSTest.(samples1, samples2) .|> pvalue
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

function AnalysisModule.get_efficiency(ks::KS, CL::Real, sampleSize::Real) 
	pValsId = isempty(findfirst( x -> x == sampleSize, ks.sampleSizes )) ? error("invalid sample size, $sampleSize, not in provided sample sizes $ks.sampleSizes") : findfirst( x -> x == sampleSize, chi2.sampleSizes ) 

	return AnalysisModule.get_efficiency(ks.pVals[pValsId], ks.CL) 
end

function get_best_sample_size(ks::KS, sampleSizes, CL, nSamples, nInARow)
	effs = get_efficiency.(ks.pVals, CL, nSamples) 
	return get_best_sample_size(effs, sampleSizes, nInARow)
end