module AnalysisModule
using FHist 
using RecipesBase, StatsBase, UnROOT, HypothesisTests

include("./Misc.jl")
export  get_min_value,			# returns the minimal value in a matrix that is !!not NaN!!
		get_min_idx,			# returns cartesian coordinate index of the minimal value of the matrix (obtained by get_min_value)
		cartesianIdx_to_range,	# converts idx to range of angles
		get_needed_statistics,	# returns Stats (as defined in PlutAnalysis.jl)
		prop_error,				# returns a confidence interval for proportion (of 100 events)
		fill_from_root_file,  	# fill an array with values from root file
		get_samples,
		get_efficiency,
		get_residuals,
		get_residuals_errors,
		get_kappa,
		DeltaM,
		Mmin 

include("./BBB.jl")
export  get_r, 					# returns the ratio of two bin heights (normalized)
		get_delta_r,			# returns the uncertainty on r (as defined in PlutAnalysis.jl)
		get_M,					# returns \tilde M, the minimal number of counts needed in the respective bin to distinguish two spectra (as defined in PlutAnalysis.jl)
		get_r_min,				# returns best case ratio
		get_r_max,				# returns worst case ratio
		get_needed_statistics,	# returns Stats (as defined in PlutoAnalysis.jl)
		get_r_map,				# returns a map of r ratios 
		get_r_max_map,			# returns a map of rmax ratios
		get_r_min_map,			# returns a map of rmin ratios
		get_delta_r_map,        # returns a map of delta r ratios
		get_min_stats_map 		# returns a map of min stats needed

export BBB 						# main Module object

include("./Chi2.jl")
export ChisqTest, get_best_sample_size, get_efficiency
export Chi2 					# main Module object

include("./KS.jl")
export get_best_sample_size, get_efficiency
export KS						# main Module object

end #MODULE END
