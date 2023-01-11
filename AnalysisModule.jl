module AnalysisModule
using FHist 

export  get_r, 					# returns the ratio of two bin heights (normalized)
		get_delta_r,			# returns the uncertainty on r (as defined in PlutAnalysis.jl)
		get_M,					# returns \tilde M, the minimal number of counts needed in the respective bin to distinguish two spectra (as defined in PlutAnalysis.jl)
		get_r_min,				# returns best case ratio
		get_r_max,				# returns worst case ratio
		get_min_value,			# returns the minimal value in a matrix that is !!not NaN!!
		get_min_idx,			# returns cartesian coordinate index of the minimal value of the matrix (obtained by get_min_value)
		cartesianIdx_to_range,	# converts idx to range of angles
		get_needed_statistics	# returns Stats (as defined in PlutAnalysis.jl)



function get_r( h1, h2, xMin, xMax, Δϕ )
	    M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) # get the sum of bincounts in the range xMin (included) - xMax (excluded)
	    N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))   
	    
	    Mnormed = M / integral(h1)
	    Nnormed = N / integral(h2)
	    
	    if( Mnormed / Nnormed <= 1 )
	        return Mnormed/Nnormed
	    else
	        return Nnormed/Mnormed
	    end
end

function get_delta_r( h1, h2, xMin, xMax, Δϕ )
	M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) # get the sum of bincounts in the range xMin (included) - xMax (excluded)
	N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))    
	
	return get_r(h1, h2, xMin, xMax, Δϕ)*sqrt(1/M + 1/N)
end


function get_M( h1, h2, xMin, xMax, Δϕ, nSigma )
	r = get_r(h1, h2, xMin, xMax, Δϕ)  
	numerator = r + sqrt(r)
	denominator = 1-r
	return nSigma^2*( numerator / denominator )^2
end

function get_M(r::Real, nSigma::Real)
	numerator = r + sqrt(r)
	denominator = 1-r
	
	return nSigma^2*( numerator / denominator )^2
end


function get_r_min( h1, h2, xMin, xMax, Δϕ )
	M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) # get the sum of bincounts in the range xMin (included) - xMax (excluded)
	N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))   
	
	Mtot = integral(h1)
	Ntot = integral(h2)
	
	εM = M / Mtot
	εN = N / Ntot
	
	ΔεM = sqrt(M) / Mtot
	ΔεN = sqrt(N) / Ntot
	
	if( εM / εN <= 1 )
		return (εM - ΔεM)/(εN + ΔεN)
	else
		return (εN - ΔεN)/(εM + ΔεM)
	end
end

function get_r_max( h1, h2, xMin, xMax, Δϕ )
	M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) # get the sum of bincounts in the range xMin (included) - xMax (excluded)
	N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))   
	
	Mtot = integral(h1)
	Ntot = integral(h2)
	
	εM = M / Mtot
	εN = N / Ntot
	
	ΔεM = sqrt(M) / Mtot
	ΔεN = sqrt(N) / Ntot
	
	if( εM / εN <= 1 )
		return (εM + ΔεM)/(εN - ΔεN)
	else
		return (εN + ΔεN)/(εM - ΔεM)
	end
end

function get_min_value(mat)
	minVal = Inf
	for i in CartesianIndices(mat)
		if( mat[i] != NaN && mat[i] < minVal  )
			minVal = mat[i]
		end
	end
	return minVal
end

function get_min_idx(mat)
	minVal = Inf
	idx = (0,0)
	for i in CartesianIndices(mat)
		if( mat[i] != NaN && mat[i] < minVal  )
			minVal = mat[i]
			idx = i
		end
	end
	return idx
end

function cartesianIdx_to_range(idx, Δϕ)
	return idx[2]*Δϕ-Δϕ, idx[1]*Δϕ
	
end

function get_needed_statistics(h1, h2, xMin, xMax, Δϕ, nSigma)
	M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) 
	N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))   
	
	Mtot = integral(h1)
	Ntot = integral(h2)
	
	εM = M / Mtot
	εN = N / Ntot
	
	Mr = get_M( h1, h2, xMin, xMax, Δϕ, nSigma )
	ε = εM < εN ? εM : εN

	return Mr / ε
end

function get_needed_statistics(h1, h2, xMin, xMax, Δϕ, nSigma, Mr)
	M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) 
	N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))   
	
	Mtot = integral(h1)
	Ntot = integral(h2)
	
	εM = M / Mtot
	εN = N / Ntot
	
	ε = εM < εN ? εM : εN
	return Mr / ε
end
	

end #MODULE END
