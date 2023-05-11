mutable struct BBB{T<:Real}
	vector1::Vector{T}
	vector2::Vector{T}
	xMin::T # minimum value of the array (i.e. for angles it would be 0 degrees)
	xMax::T # maximum value of the array (i.e. for angles it would be 180 degrees)
	stepSize::T
	nSigma::T
	minEvents::T
	minROI::Tuple{T,T}

	function BBB(vector1::Vector{T}, vector2::Vector{T}, xMin::T, xMax::T, stepSize::T, nSigma::T) where T
		if( (xMax-xMin)/stepSize%1 != 0.0 )  						# check the boundaries of the histogram
			error("Range must be integer divisible by stepsize! ")
		end
		matStatsMr::Matrix{T} 	    = get_min_stats_map(vector1, vector2, xMin, xMax, stepSize, nSigma)   # get matrix of min stats
		replace!( matStatsMr, 0.0 => NaN  ) 
		replace!( matStatsMr, Inf => NaN  ) 


		minEvents::T 			 = round(get_min_value(matStatsMr), digits = 2)                         # find minimum value
		minROIidx::CartesianIndex{2} = get_min_idx(matStatsMr)                                              # find where min value is at
		minROI::Tuple{<:T,<:T} = cartesianIdx_to_range(minROIidx, stepSize)                           # convert from index to cartesian coordinates

		new{T}(vector1::Vector{T}, vector2::Vector{T}, xMin::T, xMax::T, stepSize::T, nSigma::T, minEvents::T, minROI::Tuple{<:T,<:T}) 
	end

end

function get_r( h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, Δϕ::Real )
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

function get_delta_r( h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, Δϕ::Real  )
	M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) # get the sum of bincounts in the range xMin (included) - xMax (excluded)
	N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))    
	
	return get_r(h1, h2, xMin, xMax, Δϕ)*sqrt(1/M + 1/N)
end


function get_M( h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, Δϕ::Real, nSigma::Real )
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


function get_r_min( h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, Δϕ::Real  )
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

function get_r_max( h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, Δϕ::Real  )
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



function get_needed_statistics(h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, Δϕ::Real, nSigma::Real)
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

function get_needed_statistics(h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, Δϕ::Real, nSigma::Real, Mr::Matrix{<:Real})
	M = sum(lookup.(h1, xMin:Δϕ:xMax-Δϕ)) 
	N = sum(lookup.(h2, xMin:Δϕ:xMax-Δϕ))   
	
	Mtot = integral(h1)
	Ntot = integral(h2)
	
	εM = M / Mtot
	εN = N / Ntot
	
	ε = εM < εN ? εM : εN
	return Mr / ε
end


function get_r_map(h1::Hist1D, h2::Hist1D, xMin::Real, xMax::Real, stepSize::Real)
    matRatios = zeros(Int(xMax/stepSize),Int(xMax/stepSize))

    for minROI in xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
        for maxROI in minROI+stepSize:stepSize:xMax
            r = get_r(h1, h2, minROI, maxROI, stepSize)
    
            matRatios[Int(maxROI/stepSize), Int(minROI/stepSize)+1] = r
        end
    end

    return matRatios
end

function get_r_map(vec1::Vector{<:Real}, vec2::Vector{<:Real}, xMin::Real, xMax::Real, stepSize::Real)
    h1 = Hist1D( vec1, xMin:stepSize:xMax )
    h2 = Hist1D( vec2, xMin:stepSize:xMax )

    matRatios = zeros(Int(xMax/stepSize),Int(xMax/stepSize))

    for minROI in xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
        for maxROI in minROI+stepSize:stepSize:xMax
            r = get_r(h1, h2, minROI, maxROI, stepSize)
    
            matRatios[Int(maxROI/stepSize), Int(minROI/stepSize)+1] = r
        end
    end

    return matRatios
end

function get_r_map(bbb::BBB)
    get_r_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)
end

function get_r_max_map(vec1::Vector{<:Real}, vec2::Vector{<:Real}, xMin::Real, xMax::Real, stepSize::Real)
    h1 = Hist1D( vec1, xMin:stepSize:xMax )
    h2 = Hist1D( vec2, xMin:stepSize:xMax )

    matRatios = zeros(Int(xMax/stepSize),Int(xMax/stepSize))

    for minROI in xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
        for maxROI in minROI+stepSize:stepSize:xMax
            rmax =  get_r_max(h1, h2, minROI, maxROI, stepSize) <= 1.0 ? 
			        get_r_max(h1, h2, minROI, maxROI, stepSize) : NaN   
    
            matRatios[Int(maxROI/stepSize), Int(minROI/stepSize)+1] = rmax
        end
    end

    return matRatios
end

function get_r_max_map(bbb::BBB)
    get_r_max_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)
end

function get_r_min_map(vec1::Vector{<:Real}, vec2::Vector{<:Real}, xMin::Real, xMax::Real, stepSize::Real)
    h1 = Hist1D( vec1, xMin:stepSize:xMax )
    h2 = Hist1D( vec2, xMin:stepSize:xMax )

    matRatios = zeros(Int(xMax/stepSize),Int(xMax/stepSize))

    for minROI in xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
        for maxROI in minROI+stepSize:stepSize:xMax
            rmin =  get_r_min(h1, h2, minROI, maxROI, stepSize) 
    
            matRatios[Int(maxROI/stepSize), Int(minROI/stepSize)+1] = rmin
        end
    end

    return matRatios
end

function get_r_min_map(bbb::BBB)
    get_r_min_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)
end

function get_delta_r_map(vec1::Vector{<:Real}, vec2::Vector{<:Real}, xMin::Real, xMax::Real, stepSize::Real)
    h1 = Hist1D( vec1, xMin:stepSize:xMax )
    h2 = Hist1D( vec2, xMin:stepSize:xMax )

    matRatios = zeros(Int(xMax/stepSize),Int(xMax/stepSize))

    for minROI in xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
        for maxROI in minROI+stepSize:stepSize:xMax
            rmin = get_delta_r(h1, h2, minROI, maxROI, stepSize) 
    
            matRatios[Int(maxROI/stepSize), Int(minROI/stepSize)+1] = rmin
        end
    end

    return matRatios
end

function get_delta_r_map(bbb::BBB)
	get_delta_r_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)
end


function get_min_stats_map(vec1::Vector{<:Real}, vec2::Vector{<:Real}, xMin::Real, xMax::Real, stepSize::Real)
    h1 = Hist1D( vec1, xMin:stepSize:xMax )
    h2 = Hist1D( vec2, xMin:stepSize:xMax )

    matRatios = zeros(Int(xMax/stepSize),Int(xMax/stepSize))

    for minROI in xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
        for maxROI in minROI+stepSize:stepSize:xMax
            rmax = get_r_max(h1, h2, minROI, maxROI, stepSize) <= 1.0 ? 
			       get_r_max(h1, h2, minROI, maxROI, stepSize) : NaN  
            
            matRatios[Int(maxROI/stepSize), Int(minROI/stepSize)+1] = rmax <= 1.0 ? get_needed_statistics(h1, h2, minROI, maxROI, stepSize, nSigma) : NaN
        end
    end

    return matRatios
end

function get_min_stats_map(bbb::BBB)
    get_min_stats_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize, bbb.nSigma)
end



	