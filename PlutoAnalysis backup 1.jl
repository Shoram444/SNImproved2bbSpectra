### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 411eeb6b-6902-4cd9-8bda-8ebae1de81b1
# loading necessary packages
begin
	using Revise 
	using StatsPlots, UnROOT, StatsBase, Polynomials, ColorSchemes, Suppressor  , HypothesisTests, LaTeXStrings
	using FHist, DataFramesMeta, Distributions, DataFrames , RecipesBase
	using PlutoUI

	AM = include("./src/AnalysisModule.jl")
end

# ╔═╡ 5565017c-534f-46a5-b86b-37753c42da24
# setting plotting theme
begin 
	gr() 
	default(fmt = :jpg)
	theme(
	    :dao;
	    size           = (800, 800),
	    legend         = :topleft,
	    guidefontsize  = 16,
	    tickfontsize   = 12,
	    titlefontsize  = 16,
	    legendfontsize = 12,
	    left_margin    = 4Plots.mm,
	    right_margin   = 8Plots.mm,
	    top_margin     = 4Plots.mm,
	    bottom_margin  = 6Plots.mm,
	    dpi            = 200,
	    :colorbar_titlefontsize => 20,
	    widen = :false,
		:markerstrokewidth => 1,
		:markerstrokecolor => :black,
	);

end

# ╔═╡ afeccb00-7abf-11ed-194a-558eaaa68040
md"""
Analysis notebook for improved ``2\nu\beta\beta`` spectra. 

The **goal** of this analysis is to compare 2 (slighty) different angular distributions. 

The **methodology** used in the analysis is as follows. First, the angular distribution of ``2\nu\beta\beta`` is given by the equation:

$\frac{d\Gamma}{dcos(\theta)} \sim N(1 + K (\xi_{31}, \xi_{51})cos\theta)$,

equation (24) from paper by [Ovidiu](https://www.mdpi.com/2218-1997/7/5/147). In order to properly sample angular distributions, the parameter $K^{2\nu}$ must be provided. This parameter is defined (in eq. 26 from Ovidiu's paper) as:

$K^{2\nu} = - \frac{H_0 + \xi_{31}H_2 + \frac{5}{9}\xi_{31}^2H_{22}+ (\frac{2}{9}\xi_{31}^2 + \xi_{51})H_{4}}{G_0 + \xi_{31}G_2 + \frac{1}{3}\xi_{31}^2G_{22} + (\frac{1}{3}\xi_{31}^2 + \xi_{51})G_4}$

The single-electron energy distribution is in turn given by eq. 35 from [Fedor's publication](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.97.034315):

$\frac{d\Gamma}{dT_{e}} \sim \frac{dG_0}{dT_{e}} + \xi_{31}\frac{dG_2}{dT_{e}} + \frac{1}{3}(\xi_{31})^2\frac{dG_{22}}{dT_{e}} + (\frac{1}{3}(\xi_{31})^2 + \xi_{51})\frac{dG_4}{dT_{e}}$. 

Figure below shows ``K(\xi_{31}, \xi_{51} = 0.1397)``:


"""

# ╔═╡ 716b3e40-73e6-4d80-8126-8ec10d0d64fb
plot( 
	-1:0.001:1.0, 
	AM.get_kappa.(-1:0.001:1.0), 
	xlabel = "ξ31", 
	ylabel = "K2ν", 
	title = "fixed ξ51 = 0.1397", 
	label = "",
	lw = 4
)

# ╔═╡ 6e1c4c5b-08a3-44fa-bd34-770fd38c03a3
md"""

In the initial stage of the analysis 2 sets of ``2\nu\beta\beta`` spectra were simulated in Falaise, ``1e8`` simulated events each. First, a spectrum (*standard spectrum - SM*) with ``K = -0.88`` (taken from Decay0, which is currently in Falaise) and second spectrum is the refined one with ``Κ = -0.6639; ξ31 = 0.37``. The simulated data were passed through flreconstruct and the following data cuts were applied: 

1. 	Two negatively charged tracks reconstructed,
2. 	Two vertices on the source foil, within given distance from each other,
3. 	Sum of electron energies within the range: ``E_{sum} \in (0, 3500)~keV``,
4. 	Two individual Optical Module hits,
5. 	Two associated Optical Module hits.

From this, angular and single-electron energy distributions were obtained. The angle ``\phi`` represents the **escape** angle - ie. the angle which can be measured with SuperNEMO, the angle between the two electrons at the moment they escape the source foil. 

The two spectra are shown below. 
"""

# ╔═╡ a7e6341c-6f8b-42a6-b440-368fea258855
begin
	fName1 = "/home/shoram/Work/PhD_Thesis/Job17/Data/2vbb_Se82_Falaise_EneThetaPhi_1e8E.root"                           # root file for reference spectra
	fName2 = "/home/shoram/Work/PhD_Thesis/Job17/Data/G0-G4_xi31_037_kappa_m06639_1e8E.root"   # root file for compared spectra
	
	file1 = ROOTFile( fName1 )
	file2 = ROOTFile( fName2 )

	singleElectronEnergies1 = 
		AM.fill_from_root_file(file1, "tree","reconstructedEnergy1") .+ AM.fill_from_root_file(file1, "tree", "reconstructedEnergy2") # vector of single-electron energies for reference spectrum    
	singleElectronEnergies2 = 
		AM.fill_from_root_file(file2, "tree","reconstructedEnergy1") .+ AM.fill_from_root_file(file2, "tree", "reconstructedEnergy2") # vector of single-electron energies for compared spectrum    

	phi1 = AM.fill_from_root_file(file1, "tree","phi") # vector phi angles for reference spectrum    
	phi2 = AM.fill_from_root_file(file2, "tree","phi") # vector phi angles for compared spectrum   
end

# ╔═╡ 6c4cc9cc-3344-4ff0-bed3-df7ba1d00bc7
begin
	Δϕ = 5 # setting bin width
	ΔE = 60
end

# ╔═╡ da3a265b-6230-4414-b4c4-01abc5e9b3e8
begin
	shPhi = histogram(
		[phi1, phi2], 
		normed=:true, 
		nbins=0:Δϕ:180, 
		label=["Κ=-0.88 (SM)" "Κ=-0.6639 (refined)"], 
		xlabel="escape angle φ [°]" ,
		ylabel="normalized count rate",
		lw=4
	)

	resPhi = scatter(
		midpoints(0:Δϕ:180),
		AM.get_residuals( phi1, phi2, 0:Δϕ:180, true ),
		yerr=AM.get_residuals_errors( phi1, phi2, 0:Δϕ:180, true ),
		mc=:black,
		label="",
		ylabel="refined/SM",
		xlabel="escape angle φ [°]",
		ms=1.5,
		ylims=(0.8, 1.2)
	)

	shEne = histogram(
		[singleElectronEnergies1, singleElectronEnergies2], 
		normed=:true, 
		nbins=0:ΔE:3000, 
		label=["SM" "ξ31=0.37 (refined)"], 
		xlabel=L"\textrm{single ~electron ~energy ~T_e ~[keV]}" ,
		ylabel="normalized count rate",
		legend=:best,
		lw=4
	)

	resEne = scatter(
		midpoints(0:ΔE:3000),
		AM.get_residuals( 
			singleElectronEnergies1, singleElectronEnergies2, 0:ΔE:3000, true ),
		yerr=AM.get_residuals_errors( 
			singleElectronEnergies1, singleElectronEnergies2, 0:ΔE:3000, true ),
		mc=:black,
		label="",
		ylabel="refined/SM",
		xlabel=L"\textrm{single ~electron ~energy ~T_e ~[keV]}",
		ylims=(0.8, 1.2),
		ms=1.5
	)

	plot(
		shPhi,
		shEne, 
		resPhi,
		resEne,
		layout = grid(2,2), 
		plot_title = "reconstructed distributions",
		size = (1400, 800),
		thickness_scaling = 1.4
	)
end

# ╔═╡ 95a6b54a-a692-4477-954f-849f015e6825
md"""

Looking at the two figures for the residuals, we can notice the following. First, for both distributions ``\phi`` and Energy the biggest difference is at the low/high edges (though, for angles, we must take into account that the detector performance in such regions is lowered compared).
Second, while the difference for the ``\varphi`` distribution appears to be more pronounced than that for energy. 

"""

# ╔═╡ 876fd917-f94a-4404-b274-f376e4b41a44
md"""
In the next section three methods for spectral comparison will be presented:
1. Bin-by-bin (BBB)
2. ``\chi^2`` hypothesis test
3. Kolmogorov-Smirnov hypothesis test

The goal is to quantify the difference between the two angular distributions. We wish to provide the following answers: 
1. 	How many events (how large a statistics) must be measured to be able to distinguish two spectral shapes within given ``n_{\sigma}`` confidence?
2. 	Is it feasible to obtain such statistics within the 5 year data-taking period of SuperNEMO? If not, what would have to be the parameters?
2. 	What is the sensitivity of SuperNEMO toward the given distribution?
"""

# ╔═╡ 85cb8fd5-70c9-423a-8ba0-94a8632d5291
md"""
**Bin-By-Bin method**

The main idea arises from comparing the difference between the two distributions in terms of bin heights (``h_{i}^{j}``; ``k \in (1,2)``). First, we introduce the ratio ``r_{i} \equiv \frac{h_{i}^{j,k}}{h_{i}^{k,j}}``, where whether ``i`` is in the numerator or ``j`` is determined dependent on whether ``j < k`` or ``j > k``, respectively. Furthermore, to be more *fair* in comparisons, we exchange ``h_{i}^{j,k}`` with ``\varepsilon_{i}^{j,k} \equiv \frac{h_{i}^{j,k}}{total\_number\_of\_events}``, ie. the normalized bin height. Thus, the ratio now becomes: ``r_{i} \equiv \frac{\varepsilon_{i}^{j,k}}{\varepsilon_{i}^{k,j}}``

To determine ROI where the ratio ``r_i`` is most favourable for our purposes, that is the lowest number, we create *maps* of various ranges of ``\phi``.
The maps are created by taking some range ``\phi \in (\phi_{min}, \phi_{max})``, and calculating the respective ``r_i``. The results are shown in the figure below. Each square in the figure represents a certain range, to be read out by the upper-left corner of the square. (That is, the range ``\phi \in (0, \Delta\phi)`` is represented by the square in the down-left corner, the very first square.)

"""


# ╔═╡ 59d2eae1-3873-41be-90eb-8d989c27754d
md"""
The method is enveloped in the `BBB.jl` file as part of the `AnalysisModule.jl`. The main object for this method is the `BBB` type, which takes in the following inputs:
1. `vector1::{Vector{<:Real}}`: First (reference) vector of values, i.e. measured ``\varphi's`` from standard distribution
2. `vector2::{Vector{<:Real}}` Second (compared) vector of values, i.e. measured ``\varphi's`` from the refined distribution
3. `xMin::Real` minimum value for the binning, i.e. `0°`
4. `xMax::Real` maximum value for the binning, i.e. `180°`
5. `stepSize::Real` bin step, i.e. `5°`
6. `nSigma::Real` number of sigmas

At construction `BBB` has additional fiels: 
1. `minEvents::Real` minimum required number of events 
2. `minROI::Tuple{<:Real,<:Real}` best region of interest
"""



# ╔═╡ f3af938f-7e1d-4d7b-acda-dbf4682d5607
md"""
We first look at the behaviour of `M(r)`. Notice that for ``r \rightarrow 1.0, ~M(r)`` grows rapidly. 
"""

# ╔═╡ 8b627eb8-8d8b-431a-9768-0764045db3c8
plot(0.8:0.001:0.999, AM.Mmin.(0.8:0.001:0.999, 1e5, 1e5), 
	xlabel = "r", 
	ylabel = "M(r)",
	label = "",
	lw = 4,
	size = (800,400),
	yscale= :log10,
	yticks= [10^3, 10^4, 10^5, 10^6, 10^7, 10^8],
	xticks= 0.8:0.02:1.0,
	thickness_scaling = 1.2
)

# ╔═╡ 6c8d381d-994a-45b3-b59e-73d3c3f4e29d
begin
	BBBPhi = AM.BBB(phi1, phi2, 0, 180, Δϕ, 1)
	BBBEne = AM.BBB(singleElectronEnergies1, singleElectronEnergies2, 0, 3000, ΔE, 1)
end

# ╔═╡ f1f65c2b-480e-463b-aa3d-ba551a097843
md"""
We can look at what the r-ratios look like for each distribution
"""

# ╔═╡ af03769b-34d1-420a-8397-709b8f6218a3
begin
	rBBBPhi = AM.get_r_map(BBBPhi) 
	replace!(rBBBPhi, 0.0 => NaN)

	hmrBBBPhi = heatmap(
		midpoints(0:Δϕ:180),
		midpoints(0:Δϕ:180),
		rBBBPhi,
		c = :jet,
		xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "r",
	    right_margin = 5Plots.mm,
	    left_margin = 5Plots.mm,
	    size = (800,600),
		title="r-ratios for BBBPhi"
	)

	rBBBEne = AM.get_r_map(BBBEne) 
	replace!(rBBBEne, 0.0 => NaN)

	hmrBBBEne = heatmap(
		midpoints(0:ΔE:3000),
		midpoints(0:ΔE:3000),
		rBBBEne,
		c = :jet,
		xlabel = "Emin",
	    ylabel = "Emax",
	    colorbar_title= "r",
	    right_margin = 5Plots.mm,
	    left_margin = 5Plots.mm,
	    size = (800,600),
		title="r-ratios for BBBEne"
	)

	plot(
		hmrBBBPhi, 
		hmrBBBEne, 
		aspect_ratio = 1, 
		size = (1600,750), 
		thickness_scaling = 1.4
	)
end

# ╔═╡ ef2bb53e-c19c-490d-a116-c3cba8a84c25
md"""
This map should be read in the following way:
Each square represents a specific ROI defined by (min, max) value to which the square corresponds. (i.e. the very first square - left, bottom corner - of φ distribution corresponds to ROI: φ ∈ (0,5)°. The one above it then corresponds to ROI: φ ∈ (0,10)°, and so on. )
"""

# ╔═╡ a5037ddd-09bd-4f59-b9b8-c23dba5de9b8
md"""
In the figure we can see the calculated ``r_i`` for each range of ``\phi``'s and Energies. The most ideal value is the lowest - greatest difference between the two spectra.

This however, does not tell the whole story yet. We must consider the uncertainty of the result as well as the uncertainty of the simulation. Furthermore, since the aim was to answer **how many events are required** to distinguish the two spectra, we still have some steps to take. 

We will therefore produce a few more calculations. 
"""

# ╔═╡ 31f648a2-cdb3-4def-8f04-dbfc65787536
md"""
First, to calculate the number of events required, we begin with the following:

Assume ``M`` represents the number of events in the *smaller* bin i and ``N`` represents the number of events in the *larger* bin. Then in order to distinguish the two bins from each other at the level of ``n_{\sigma}``, their difference must be at least equal (or greater) than the sum of their respective uncertainties. 

``M - N = n_{\sigma}( \Delta M + \Delta N ); \Delta M = \sqrt M, \Delta N = \sqrt N``

Dividing the equation by ``M``, substituting ``r = \frac{M}{N}`` and rearranging yields:

``\tilde M(r) = n_{\sigma}^2(\frac{r+\sqrt{r}}{1-r})^2``,
where ``\tilde M`` simply represents the minimum number of events in ``M`` required to distinguish ``M`` from ``N`` by ``n_{\sigma}``. 

Now, using uncertainty propagation, we can find the uncertainty on ``\tilde M`` as:

``\Delta \tilde M = n_{\sigma}^2(\frac{r+\sqrt{r}}{1-r})^3r\sqrt{1/M + 1/N}``. 

**The uncertainty in ``\tilde M`` can be improved by obtaining higher statistics**. 

**However, due to the fact that we do not precisely understand the detector angular correlations, we do not know precisely the analytical value for ``r`` either**. (If we knew perfectly the correlation ``\theta \rightarrow \phi``, we could obtain ``r`` analytically from the input angular distributions). We must, therefore, take into account the uncertainty on ``r``.  

``\Delta r = r\sqrt{1/M + 1/N}; ``

Since the behavior of r is very much non-linear (``\tilde M`` rapidly explodes for ``r`` close to 1.0) we consider the maximum (and minimum) uncertainties:

``r_{max} = \frac{M + \Delta M }{N - \Delta N}``
``r_{min} = \frac{M - \Delta M}{N - \Delta N}``.

Here, the worst case scenario (``r_{max}``) can in some cases exceed 1.0, (ie. ``r > 1.0``). Such regions will be removed from the analysis and will show in the maps as empty squares. 

To answer the question *how many events SN needs to measure to distinguish* we can convert the ``\tilde M(r)`` into the number of events needed ``S`` as follows:
``S = \tilde M(r) / \varepsilon ``. That is, we scale the number of needed events in the bin by the proportion of the total events that bin represents. ``S`` then gives the total statistics needed to obtain desired ``\tilde M(r)``.

Finally, let us look at all of the mentioned values. 
In the figure below we show maps for ``S (r)``. 
"""

# ╔═╡ 21771deb-b4a7-4b08-aba6-46f2f7089788
begin
	sBBBPhi = AM.get_min_stats_map(BBBPhi) 
	replace!(sBBBPhi, 0.0 => NaN)

	hmsBBBPhi = heatmap(
		midpoints(0:Δϕ:180),
		midpoints(0:Δϕ:180),
		sBBBPhi,
		c = :jet,
		xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "\nStats needed",
	    right_margin = 5Plots.mm,
	    left_margin = 5Plots.mm,
	    size = (800,600),
		title="Stats needed for BBBPhi",
		colorbar_scale = :log10
	)

	sBBBEne = AM.get_min_stats_map(BBBEne) 
	replace!(sBBBEne, 0.0 => NaN)

	hmsBBBEne = heatmap(
		midpoints(0:ΔE:3000),
		midpoints(0:ΔE:3000),
		sBBBEne,
		c = :jet,
		xlabel = "Emin",
	    ylabel = "Emax",
	    colorbar_title= "\nStats needed",
	    right_margin = 15Plots.mm,
	    left_margin = 5Plots.mm,
	    size = (800,600),
		title="Stats needed for BBBEne",
		colorbar_scale = :log10
	)

	plot(
		hmsBBBPhi, 
		hmsBBBEne, 
		size = (1600,700), 
		thickness_scaling = 1.4
	)
end

# ╔═╡ 6c3b9653-3d26-4c07-a362-de5649879d12
md"""
First of all, the *excluded* regions are the ones for which the *maximum-error* exceedes physical values (negative under square-root).

We can see in these maps, that the values for *stats-needed* can take a very wide range. The best values, for `nSigma = 1`, are found as:
1. φ distribution: $(BBBPhi.minEvents) for ROI: $(BBBPhi.minROI) degrees
2. Energy distribution: $(BBBEne.minEvents) for ROI: $(BBBEne.minROI) keV
"""

# ╔═╡ ed4d68ee-131a-4a44-99a3-eac9ece030c8
md"""
The concludions from this method should be taken with a grain of salt as it is very non-standard. We therefore use more common approaches: ``\chi^2`` and KS.

**The following text will deal with ``\chi^2`` and `KS` methods together**
"""

# ╔═╡ ac6ee28c-72a2-4e79-8cb1-2f5a60eea8c2
md"""
The goal is to find the ideal value of number of events needed to distinguish (`S`).
To do so, we construct ``H_0 \equiv ~tested~ distributions ~have ~the ~same ~underlying ~distribution``. That is, within null hypothesis we expect the two compared spectra to be the same. Thus our condition is to reject this hypothesis (we want ``p-value < \alpha``).

**The methodology is as follows**: 
1. The spectra are split into `100` random subsets of various sizes `M`
2. For each subset a KS and ChiSquare hypothesis test is performed and `p-value` is extracted. We obtain 100 p-values for each `M`. 
3. We define efficiency ``\varepsilon \equiv N_{rejected ~H_0}/N_{total = 100}`` 
3. For various values of CL (i.e. 90%, 95%) the corresponding `S` is found as the minimum `M` for which the efficiency of rejecting ``H_0`` is **100% for three consequitive times** for increasing `M's`. 
"""

# ╔═╡ 53484f4e-129e-4aa1-8e24-959eaa9cdb3d
md"""
We define `sampleSizes` (M) to be sizes from 10000 events to 19000 events (in step 1000) and from 20000 to 150000 (in steps of 10000) and up to 400000 in steps of 50000.

The methods are, again, contained within their own submodules: `Chi2.jl` and `KS.jl`. 
The input arguments for both are the same:
1. `vector1::Vector{T}`
2. `vector2::Vector{T}`
3. `xMin::Real` # minimum value of the array (i.e. for angles it would be 0 degrees)
4. `xMax::Real` # maximum value of the array (i.e. for angles it would be 180 degrees)
5. `stepSize::Real`
6. `sampleSizes::Vector{<:Real}`
7. `CL::Real`

The constructor then makes three more fieds: 
1. `pVals::Vector{Vector{<:Real}}` - container for p-values, for each subset of size M we get 100 p-values
2. `efficiencies::Vector{<:Real}` - vector of efficiencies for each subset size
3. `minEvents::Real` - minimum events required
"""

# ╔═╡ 3d268ae3-f1db-489d-887d-de555417c467
begin
	sampleSizes = vcat(collect(20_000:10_000:100_000), collect(150_000:50_000:800_000) )
	xticks=(1:length(sampleSizes), sampleSizes)
	CL = 0.95
end

# ╔═╡ e7e5a693-efe1-4a2d-8aa1-4c24ee377168
md"""
**Disclaimer: For ``\chi^2`` there is an extra condition that must be filfilled taking into account the methodology. The binning must be chosen so that at least 5 events are in each bin. If this is not fulfilled, it will throw error.**
"""

# ╔═╡ 852b33c3-99fb-4d3f-b22c-737a301a9cd8
md"""
Above is a plot of mean p-values for each combination. We can see that the `energy` distribution needs more statistics to reject ``H_0``, in fact we have not even reach enough to obtain desired statistics. We look a bit closer at the data by looking at a boxplot, which shows the median and 2nd, 3rd quantile and outliers. 
"""

# ╔═╡ 51c2801e-7051-490a-bd68-1dceb9d7f7f9
md"""
First, we can look at a boxplot of the obtained p-values obtained from KS test for ϕ for each sample size. 
We can see that for smaller sample sizes the p-values vary greatly. Whereas for larger sample sizes, i.e. ``N>100k`` events it is very close to zero.
"""

# ╔═╡ 05700042-ca64-40d7-8ef2-76a75897fe11
md"""
We can also look at the other tests: 
1. ϕ distribution using KS test
2. Single electron spectra (E) using KS test
3. ϕ distribution using Chi2 test
3. E distribution using Chi2 test

Right away we can see that KS test is more strict when it comes to rejecting H0. Furthermore, it seems that the angular distribution is the more promising candidate for distinguishing the spectra.
"""

# ╔═╡ 75f2bc28-529e-47e1-95cf-e8a6b2f864ec
md"""
However, since it can be seen that the p-values vary greatly, we can perform an additional analysis by looking at the efficiency of rejecting H0 for various sample sizes. That is, for each subset of size ``N`` we compute efficiency as ``\varepsilon = \frac{N_{reject}}{N}``. 
"""

# ╔═╡ 5163e320-0711-42e8-89c6-90994f612c56
md"""
Finally we can create a table of *S* for the various defined methods and sample sizes.
"""

# ╔═╡ 15149b75-4ce7-472c-adfa-8d4ad153ba95
md"""
CL\method | KSphi |KSEne | Χ²Phi | Χ²Ene | BBBPhi | BBB Ene
----------|-------|------|-------|-------|--------|--------
68%       |70k    |450k  |60k    |350k   |11k     |60k
90%       |150k   |750k  |60k    |350k   |30k     |160k
95%       |200k   |750k  |200k   |500k   |45k     |230k

"""

# ╔═╡ 9a2dbe13-66c9-4fa5-8404-2872601bafef

	function get_delta_r_map(vec1::Vector{<:Real}, vec2::Vector{<:Real}, xMin::Real, xMax::Real, stepSize::Real)
	    h1 = Hist1D( vec1, xMin:stepSize:xMax )
	    h2 = Hist1D( vec2, xMin:stepSize:xMax )
	
	    matRatios = zeros(Int(xMax/stepSize),Int(xMax/stepSize))
	
	    for minROI in xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
	        for maxROI in minROI+stepSize:stepSize:xMax
	            rmin =  AM.get_delta_r(h1, h2, minROI, maxROI, stepSize) 
	    
	            matRatios[Int(maxROI/stepSize), Int(minROI/stepSize)+1] = rmin
	        end
	    end
	
	    return matRatios
	end
	
	

# ╔═╡ b675e1b6-320c-4db5-b527-068e6e571c84
get_delta_r_map(bbb::AM.BBB) =  get_delta_r_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)

# ╔═╡ 45b56c40-2c2b-4725-9cf8-85958947f25b
get_r_min_map(bbb::AM.BBB) =  AM.get_r_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)

# ╔═╡ 9ee70dbc-e29a-4bfc-a5f8-84afe88eb00b
figDir = "/media/shoram/Extra SSD/CernBox/Work/Presentations/2023/2023_05_11_AnalysisMeeting/figs"

# ╔═╡ 02e4fab7-c042-4f80-a616-2f318f177302
begin 
	a = get_r_min_map(BBBPhi); replace!(a, 0.0 => NaN)
	title = "r_min"
	afigtitle = "r_min_phi.png"

	b = get_r_min_map(BBBEne); replace!(b, 0.0 => NaN)
	title = "r_min"
	bfigtitle = "r_min_ene.png"

	
	heatmap(
		midpoints(0:Δϕ:180),
		midpoints(0:Δϕ:180),
		a,
		c = :jet,
		xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    right_margin = 5Plots.mm,
	    left_margin = 5Plots.mm,
	    size = (800,600),
		title=title,
		# colorbar_scale = :log10
	)

	savefig(joinpath(figDir, afigtitle))


	
	heatmap(
		midpoints(0:ΔE:3000),
		midpoints(0:ΔE:3000),
		b,
		c = :jet,
		xlabel = "Emin",
	    ylabel = "Emax",
	    right_margin = 5Plots.mm,
	    left_margin = 5Plots.mm,
	    size = (800,600),
		title=title,
		colorbar_scale = :log10
	)

	savefig(joinpath(figDir, bfigtitle))
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
DataFramesMeta = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
FHist = "68837c9b-b678-4cd5-9925-8a54edc8f695"
HypothesisTests = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"
Suppressor = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
UnROOT = "3cd96dde-e98d-4713-81e9-a4a1b0235ce9"

[compat]
ColorSchemes = "~3.20.0"
DataFrames = "~1.4.4"
DataFramesMeta = "~0.12.0"
Distributions = "~0.25.79"
FHist = "~0.9.0"
HypothesisTests = "~0.10.11"
LaTeXStrings = "~1.3.0"
PlutoUI = "~0.7.50"
Polynomials = "~3.2.1"
RecipesBase = "~1.3.3"
Revise = "~3.5.0"
StatsBase = "~0.33.21"
StatsPlots = "~0.15.4"
Suppressor = "~0.2.1"
UnROOT = "~0.8.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "ff9bffc15194090e4a842a319e3e9e23ff75cc69"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArraysOfArrays]]
deps = ["Adapt", "ChainRulesCore", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "c59b725b0aadf7df93fb3de05b5e1b14029af2da"
uuid = "65a8f2f4-9b39-5baf-92e2-a9cc46fdf018"
version = "0.6.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BayesHistogram]]
git-tree-sha1 = "5d5dda960067751bc1534aba765f771325044501"
uuid = "000d9b38-65fe-4c81-bdb9-69f01f102479"
version = "1.0.7"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "64df3da1d2a26f4de23871cd1b6482bb68092bd5"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.3"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "0e5c14c3bb8a61b3d53b2c0620570c332c8d0663"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.2.0"

[[deps.CodecLz4]]
deps = ["Lz4_jll", "TranscodingStreams"]
git-tree-sha1 = "59fe0cb37784288d6b9f1baebddbf75457395d40"
uuid = "5ba52731-8f18-5e0d-9241-30f10d1ec561"
version = "0.4.0"

[[deps.CodecXz]]
deps = ["Libdl", "TranscodingStreams", "XZ_jll"]
git-tree-sha1 = "82c4c000edf64b6bda6766377e69a1028f3549ee"
uuid = "ba30903b-d9e8-5048-a5ec-d1f5b0d4b47b"
version = "0.7.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.CodecZstd]]
deps = ["CEnum", "TranscodingStreams", "Zstd_jll"]
git-tree-sha1 = "849470b337d0fa8449c21061de922386f32949d9"
uuid = "6b39b394-51ab-5f42-8807-6242bab2b4c2"
version = "0.7.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "00a2cccc7f098ff3b66806862d275ca3db9e6e5a"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.5.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

[[deps.DataFramesMeta]]
deps = ["Chain", "DataFrames", "MacroTools", "OrderedCollections", "Reexport"]
git-tree-sha1 = "a70c340c1306febfd770a932218561b5e19cf0f6"
uuid = "1313f7d8-7da2-5740-9ea0-a2ca25f37964"
version = "0.12.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "a7756d098cbabec6b3ac44f369f74915e8cfd70a"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.79"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FHist]]
deps = ["BayesHistogram", "LinearAlgebra", "MakieCore", "Measurements", "RecipesBase", "Requires", "Statistics", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "f65b3f4e8f215acd56af528bf2d8c1bde1dae65e"
uuid = "68837c9b-b678-4cd5-9925-8a54edc8f695"
version = "0.9.0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "9a0472ec2f5409db243160a8b030f94c380167a3"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "6872f5ec8fd1a38880f027a26739d42dcda6691f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.2"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "387d2b8b3ca57b791633f0993b31d8cb43ea3292"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.71.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "5982b5e20f97bff955e9a2343a14da96a746cd8c"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.71.3+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "fd9861adba6b9ae4b42582032d0936d456c8602d"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.6.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.HypothesisTests]]
deps = ["Combinatorics", "Distributions", "LinearAlgebra", "Random", "Rmath", "Roots", "Statistics", "StatsBase"]
git-tree-sha1 = "ae3b6964d58df11984d22644ce5546eaf20fe95d"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.10.11"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "72ab280d921e8a013a83e64709f66bc3e854b2ed"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.20"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LRUCache]]
git-tree-sha1 = "d862633ef6097461037a00a13f709a62ae4bdfdd"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.4.0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibDeflate]]
deps = ["libdeflate_jll"]
git-tree-sha1 = "4eec2ec0493906e321dc5e246820045d72539f98"
uuid = "9255714d-24a7-4b30-8ea3-d46a97f7e13b"
version = "0.4.1"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "946607f84feb96220f480e0422d3484c49c00239"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.19"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LorentzVectors]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "c479c9d097a35e331ef23df62febbaa23de8af27"
uuid = "3f54b04b-17fc-5cd4-9758-90c048d965e3"
version = "0.4.0"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "c1885d865632e7f37e5a1489a164f44c54fb80c9"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.5.2"

[[deps.MarchingCubes]]
deps = ["SnoopPrecompile", "StaticArrays"]
git-tree-sha1 = "ffc66942498a5f0d02b9e7b1b1af0f5873142cdc"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "12950d646ce04fb2e89ba5bd890205882c3592d7"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.8.0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Memoization]]
deps = ["MacroTools"]
git-tree-sha1 = "55dc27dc3d663900d1d768822528960acadc012a"
uuid = "6fafb56a-5788-4b4e-91ca-c0cea6611c73"
version = "0.1.14"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mixers]]
deps = ["MacroTools"]
git-tree-sha1 = "58ec7ac60dad6e8ca4553225251dfd380e3930dd"
uuid = "2a8e4939-dab8-5edc-8f64-72a8776f13de"
version = "0.1.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "efe9c8ecab7a6311d4b91568bd6c88897822fabe"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "f71d8950b724e9ff6110fc948dff5a329f901d64"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.8"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "df6830e37943c7aaa10023471ca47fb3065cc3c4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6e9dba33f9f2c44e08a020b0caf6903be540004"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.19+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "6466e524967496866901a78fca3f2e9ea445a559"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.2"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "5b7690dd212e026bbab1860016a6601cb077ab66"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.2"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SnoopPrecompile", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "02ecc6a3427e7edfff1cebcf66c1f93dd77760ec"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.1"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "6ea39b2399c92b83036ef26d8bab9cd017b9a8c4"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.1"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "97aa253e65b784fd13e83774cadc95b38011d734"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.6.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase", "SnoopPrecompile"]
git-tree-sha1 = "e974477be88cb5e3040009f3767611bc6357846f"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.11"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "fd5dba2f01743555d8435f7c96437b29eae81a17"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "82362f2a4f756951f21ebb3ac2aed094c46a5109"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.12"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "efd23b378ea5f2db53a55ae53d3133de4e080aa9"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.16"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "6954a456979f23d05085727adb17c4551c19ecd1"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.12"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "ab6083f09b3e617e34a956b43e9d51b824206932"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.1.1"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "DataValues", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "e0d5bc26226ab1b7648278169858adcfbd861780"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.4"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.Suppressor]]
git-tree-sha1 = "c6ed566db2fe3931292865b966d6d140b7ef32a9"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "e4bdc63f5c6d62e80eb1c0043fcc0360d5950ff7"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.10"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.UnROOT]]
deps = ["AbstractTrees", "ArraysOfArrays", "CodecLz4", "CodecXz", "CodecZstd", "HTTP", "IterTools", "LRUCache", "LibDeflate", "LorentzVectors", "Memoization", "Mixers", "Mmap", "Parameters", "PrettyTables", "SentinelArrays", "StaticArrays", "StructArrays", "Tables", "xrootdgo_jll"]
git-tree-sha1 = "1efb6662dea111f81dceab19335b739ee1d50f88"
uuid = "3cd96dde-e98d-4713-81e9-a4a1b0235ce9"
version = "0.8.21"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorSchemes", "ColorTypes", "Contour", "Crayons", "Dates", "FileIO", "FreeType", "LinearAlgebra", "MarchingCubes", "NaNMath", "Printf", "Requires", "SnoopPrecompile", "SparseArrays", "StaticArrays", "StatsBase", "Unitful"]
git-tree-sha1 = "f7d94d8025a0b231b9f04fa8e81a42e62e135084"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "3.3.2"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d670a70dd3cdbe1c1186f2f17c9a68a7ec24838c"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.12.2"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7928d348322698fb93d5c14b184fdc176c8afc82"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.2.9+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libdeflate_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0d17dd704e500cd69222fb8bd8ba63513eeaeb44"
uuid = "46979653-d7f6-5232-b59e-dd310c4598de"
version = "1.13.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"

[[deps.xrootdgo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3d51ea273d3cb8fe16fe3b80fea3f88a7c000cd0"
uuid = "9d84c17e-11f2-50ef-8cc9-e9701362097f"
version = "0.31.1+0"
"""

# ╔═╡ Cell order:
# ╠═411eeb6b-6902-4cd9-8bda-8ebae1de81b1
# ╠═5565017c-534f-46a5-b86b-37753c42da24
# ╠═afeccb00-7abf-11ed-194a-558eaaa68040
# ╟─716b3e40-73e6-4d80-8126-8ec10d0d64fb
# ╠═6e1c4c5b-08a3-44fa-bd34-770fd38c03a3
# ╠═a7e6341c-6f8b-42a6-b440-368fea258855
# ╠═6c4cc9cc-3344-4ff0-bed3-df7ba1d00bc7
# ╠═da3a265b-6230-4414-b4c4-01abc5e9b3e8
# ╟─95a6b54a-a692-4477-954f-849f015e6825
# ╟─876fd917-f94a-4404-b274-f376e4b41a44
# ╟─85cb8fd5-70c9-423a-8ba0-94a8632d5291
# ╟─59d2eae1-3873-41be-90eb-8d989c27754d
# ╟─f3af938f-7e1d-4d7b-acda-dbf4682d5607
# ╠═8b627eb8-8d8b-431a-9768-0764045db3c8
# ╠═6c8d381d-994a-45b3-b59e-73d3c3f4e29d
# ╟─f1f65c2b-480e-463b-aa3d-ba551a097843
# ╠═af03769b-34d1-420a-8397-709b8f6218a3
# ╟─ef2bb53e-c19c-490d-a116-c3cba8a84c25
# ╟─a5037ddd-09bd-4f59-b9b8-c23dba5de9b8
# ╟─31f648a2-cdb3-4def-8f04-dbfc65787536
# ╠═21771deb-b4a7-4b08-aba6-46f2f7089788
# ╟─6c3b9653-3d26-4c07-a362-de5649879d12
# ╟─ed4d68ee-131a-4a44-99a3-eac9ece030c8
# ╟─ac6ee28c-72a2-4e79-8cb1-2f5a60eea8c2
# ╟─53484f4e-129e-4aa1-8e24-959eaa9cdb3d
# ╠═3d268ae3-f1db-489d-887d-de555417c467
# ╟─e7e5a693-efe1-4a2d-8aa1-4c24ee377168
# ╠═f4702938-0f46-4448-8892-fd0998998a5a
# ╠═581e7415-2483-484c-825d-f321400bbb4c
# ╠═15d6fa2a-d3d4-4848-af09-7c081dd8f4e7
# ╟─9569f6b7-7773-4668-b19a-341ceb78e5c2
# ╟─852b33c3-99fb-4d3f-b22c-737a301a9cd8
# ╟─51c2801e-7051-490a-bd68-1dceb9d7f7f9
# ╟─f7f65fe2-f2cd-413d-8b6d-dbc5298617df
# ╟─05700042-ca64-40d7-8ef2-76a75897fe11
# ╠═51f461ac-be1f-4805-be42-9dde9a22edd6
# ╟─75f2bc28-529e-47e1-95cf-e8a6b2f864ec
# ╠═0c7d4fb7-623a-47bd-bf88-af88591b95d7
# ╠═0a7d15bb-c7e0-48dd-bd0d-043dc38de52e
# ╟─07f3f0b1-6d65-483d-bc0a-e59a03a4e55f
# ╟─5163e320-0711-42e8-89c6-90994f612c56
# ╠═c6c6a82c-2487-42bc-84e9-155346a41028
# ╠═15149b75-4ce7-472c-adfa-8d4ad153ba95
# ╠═9a2dbe13-66c9-4fa5-8404-2872601bafef
# ╠═b675e1b6-320c-4db5-b527-068e6e571c84
# ╠═45b56c40-2c2b-4725-9cf8-85958947f25b
# ╠═02e4fab7-c042-4f80-a616-2f318f177302
# ╠═9ee70dbc-e29a-4bfc-a5f8-84afe88eb00b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
