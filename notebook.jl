### A Pluto.jl notebook ###
# v0.19.17

using Markdown
using InteractiveUtils

# ╔═╡ 411eeb6b-6902-4cd9-8bda-8ebae1de81b1
# loading necessary packages
begin
	using Pkg
	Pkg.activate("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/")

	using Revise
	using StatsPlots, UnROOT, StatsBase, Polynomials  
	using FHist, DataFramesMeta, Distributions 

	ENV["COLUMNS"] = 2000
	ENV["LINES"] = 20	

	include("/home/shoram/Work/PhD_Thesis/SNAngularCorrelation/AngularCorrelations/MiscFuncs.jl")
	using .MiscFuncs
end

# ╔═╡ afeccb00-7abf-11ed-194a-558eaaa68040
md"""
Analysis notebook for improved ``2\nu\beta\beta`` spectra. 

The **goal** of this analysis is to compare 2 (slighty) different angular distributions. 

The **methodology** used in the analysis is as follows. First, the angular distribution of ``2\nu\beta\beta`` is given by the equation:
``\frac{d\Gamma}{dcos(\theta)} \sim N(1 + K (\psi_{31}, \psi_{51})cos\theta)``, equation (24) from paper by Odiviu. In the initial stage of the analysis 2 sets of ``2\nu\beta\beta`` spectra were simulated in Falaise, ``1e8`` simulated events each. First, a spectrum with ``K = -1.0`` and second with ``K = -0.65`` (a somewhat arbitrary value, close to what is mentioned in table (3) in Ovidiu's paper). The simulated data were passed through flreconstruct and the following data cuts were applied: 

1. 	Two negatively charged tracks reconstructed,
2. 	Two vertices on the source foil, within given distance from each other,
3. 	Sum of electron energies within the range: ``E_{sum} \in (0, 3500)~keV``,
4. 	Two individual Optical Module hits,
5. 	Two associated Optical Module hits.

From this two angular distributions were obtained. First, the so-called ``\theta`` angular distribution, where ``\theta`` represents the **decay** angle between the two electrons. ``\theta`` distribution cannot be experimentally measured. Second, the so-called ``\phi`` angular distribution, where ``\phi`` represents the **escape** angle - ie. the angle which can be measured with SuperNEMO, the angle between the two electrons at the moment they escape the source foil. 

The two spectra are shown below. 
"""

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
	    widen = :false
	);

end

# ╔═╡ 68624ab0-a258-4132-8e5a-54da853d856d
# setting paths to files, data, figures
begin
	baseDir = "/home/shoram/Work/PhD_Thesis/Job17/"

	dataDir = "/home/shoram/Work/PhD_Thesis/Job17/Data/kappa_-0.6_to_-0.7_0.01step/"
	figDir = baseDir*"Figures/"
end

# ╔═╡ 66154af6-bb40-46d8-84c4-3e350acb8676
dfs = DataFrame[]

# ╔═╡ a7e6341c-6f8b-42a6-b440-368fea258855
files= ["/home/shoram/Work/PhD_Thesis/Job17/Data/kappa_-0.65/MomPosEneThetaPhi_97e6E_k_m065.root", 
        "/home/shoram/Work/PhD_Thesis/Job17/Data/kappa_-1.0/MomPosEneThetaPhi_98e6E.root"]

# ╔═╡ 6424e544-8aae-4e2e-82ce-f07c0c9e8083
# loading data from files to a DataFrame
for (i,file) in enumerate(files)
    f = ROOTFile( file )

    angles = LazyTree(f, 
		"tree", 
		["theta", "phi", "reconstructedEnergy1", "reconstructedEnergy2"]
	) |> DataFrame
	
    kappa = i == 1 ? -0.65 : -1.0
    
    @transform! angles :kappa = kappa
    @transform! angles :ESum = :reconstructedEnergy1 + :reconstructedEnergy2
    
    push!(dfs, angles)
end

# ╔═╡ cad5c29a-f701-478d-902a-1a523137b5ee
# initializing two empty plots, one for ``\theta`` one for ``\phi``
begin
	hist1DPhi = plot();
	hist1DTheta = plot();
	nothing
end

# ╔═╡ 6c4cc9cc-3344-4ff0-bed3-df7ba1d00bc7
Δϕ = 5 # setting bin width

# ╔═╡ 739884af-76d1-46e8-8ee4-332c9baf40b8
for i in 1:length(dfs)
    κ = dfs[i].kappa[1]
    hphi = Hist1D(dfs[i].phi, 0:Δϕ:180)
    plot!(
		hist1DPhi, 
		hphi, 
		label = "κ = $κ", 
        ylabel = "counts ", 
		xlabel = "ϕ [°]", 
		lw = 2, 
		norm = :false, 
		xlims=(0,180), 
		alpha = 0.2,
    )
	
    htheta = Hist1D(dfs[i].theta, 0:Δϕ:180)
    plot!(
		hist1DTheta, 
		htheta, 
		label = "κ = $κ",
        ylabel = "counts ", 
		xlabel = "θ [°]", 
		lw = 2, 
		norm = :false, 
		xlims = (0,180), 
		alpha = 0.2
    )
end

# ╔═╡ 7edf09c5-bf7d-4668-a546-279ba3bc23d7
plot(
    hist1DPhi,
    hist1DTheta,
    size = (1200, 600),
    top_margin=10Plots.mm,
    left_margin = 10Plots.mm,
)   

# ╔═╡ 9afbcc00-97dc-4497-84a4-72a2a86cd16b
md"""
We can see the two angular distributions, ``\phi`` (right), ``\theta`` (left) for the two values of ``K``, with ``K = -0.65`` in *red* and ``K = -1.00`` in *blue*. 

We can see that, with ``\theta`` being the input (theoretical) spectrum, there is quite a bit of a drastic change moving toward ``\phi`` distribution. The individual features of the spectrum and hyptheses as to why the shape is as is are outlined in **Detector Effects** chapter. 

To get a better overview of the differences between the two spectra, let us look at the plots of residuals, defined as: ``\frac{h_{i}^2 - h_{i}^1}{h_{i}^1}*100``, where ``h_{i}^1`` represents the number of events in the *i-th* bin of the *reference* distribution (arbitrarily, distribution for both ``\phi`` and ``\theta`` was chosen to be the one with ``K = -0.65``), and ``h_{i}^2`` is the number of events in the *compared* distribution. 
"""

# ╔═╡ 5113de94-8b41-4f86-99fd-4382142edc28
# initializing reference histograms
begin
	refBinContentPhi = bincounts(Hist1D(dfs[1].phi, 0:Δϕ:180))
	refBinContentTheta = bincounts(Hist1D(dfs[1].theta, 0:Δϕ:180))
end

# ╔═╡ f06ff5ff-6890-482a-b461-bcb2d3c23fb5
# initializing empty plots for the residuals
begin
	residualsPhi = plot();
	residualsTheta = plot();
	nothing 
end

# ╔═╡ 26eb6233-2fb5-471a-bcf7-bf6e014fbc0c
for i in 1:length(dfs)
    κ = dfs[i].kappa[1]
    hPhi = Hist1D(dfs[i].phi, 0:Δϕ:180)
    binContentPhi = bincounts(hPhi)
    binCentersPhi = bincenters(hPhi)
    
    hTheta = Hist1D(dfs[i].theta, 0:Δϕ:180)
    binContentTheta = bincounts(hTheta)
    binCentersTheta = bincenters(hTheta)

	plot!(
        residualsPhi, 
        binCentersPhi,
        ( binContentPhi .- refBinContentPhi) ./ refBinContentPhi .* 100, 
        xlabel = "ϕ",
        ylabel = " (ϕi - ϕ1)/ϕ1*100 [%]",
        label  = "",
        legend = :top,
        lw= 2,
        seriestype =:stepmid
    )
    
    plot!(
        residualsTheta, 
        binCentersTheta,
        ( binContentTheta .- refBinContentTheta) ./ refBinContentTheta .* 100, 
        xlabel = "θ",
        ylabel = " (θi - θ1)/θ1*100 [%]",
        label  = "",
        legend = :bottom,
        lw= 2,
        seriestype =:stepmid
    )
end	

# ╔═╡ e2595477-a045-4330-9ef9-a246be1a3112
plot( 
    hist1DPhi,
    hist1DTheta,
    residualsPhi,
    residualsTheta,
    right_margin = 12Plots.mm,
    left_margin = 12Plots.mm,
    size = (1600,800),
    layout = @layout [ a b ; c{0.3h} d ]
)

# ╔═╡ 95a6b54a-a692-4477-954f-849f015e6825
md"""
**Disclaimer: for the moment we are dealing with unnormalized distributions (though the numbers are similar so it is not a big change after normalization)**

Looking at the two figures for the residuals, we can notice the following. First, for both distributions ``\phi`` and ``\theta`` the biggest difference is at low angles (thougg, we must take into account that the detector performance at such angles is lowered compared to higher angles).
Second, while the difference for the ``\theta`` distribution can be of the order of ``\sim 75 \%``, the change is nowhere near as drastic for the ``\phi`` distribution. This provides a challenge, as ``\theta`` is not experimentally measurable so we must deal with the smaller difference of the two. 

**From now on, we will only deal with ``\phi`` distribution**
"""

# ╔═╡ 85cb8fd5-70c9-423a-8ba0-94a8632d5291
md"""
**Methodology for statistical analysis**

The goal is to quantify the difference between the two angular distributions. We wish to provide the following answers: 
1. 	How many events (how large a statistics) must be measured to be able to distinguish two spectral shapes within given ``n_{\sigma}`` confidence?
2. 	Is it feasible to obtain such statistics within the 5 year data-taking period of SuperNEMO? If not, what would have to be the parameters?
2. 	What is the sensitivity of SuperNEMO toward the given distribution?

To fulfill the first goal, we have devised a methodology with the aim to find ROI, where the difference is statistically most significant. 

The main idea arises from comparing the difference between the two distributions in terms of bin heights (``h_{i}^{j}``; ``k \in (1,2)``). First, we introduce the ratio ``r_{i} \equiv \frac{h_{i}^{j,k}}{h_{i}^{k,j}}``, where whether ``i`` is in the numerator or ``j`` is determined dependent on whether ``j < k`` or ``j > k``, respectively. Furthermore, to be more *fair* in comparisons, we exchange ``h_{i}^{j,k}`` with ``\varepsilon_{i}^{j,k} \equiv \frac{h_{i}^{j,k}}{total\_number\_of\_events}``, ie. the normalized bin height. Thus, the ratio now becomes: ``r_{i} \equiv \frac{\varepsilon_{i}^{j,k}}{\varepsilon_{i}^{k,j}}``

To determine ROI where the ratio ``r_i`` is most favourable for our purposes, that is the lowest number, we create *maps* of various ranges of ``\phi``.
The maps are created by taking some range ``\phi \in (\phi_{min}, \phi_{max})``, and calculating the respective ``r_i``. The results are shown in the figure below. Each square in the figure represents a certain range, to be read out by the upper-left corner of the square. (That is, the range ``\phi \in (0, \Delta\phi)`` is represented by the square in the down-left corner, the very first square.)

"""

# ╔═╡ b7d5abe0-a806-4a6d-acf5-c0f3b9a75d33
# Declaring function that are used in the analysis. Their function is explained in the text.

begin
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
	
	
	function get_M_min( h1, h2, xMin, xMax, Δϕ, nSigma )
	    r = get_r(h1, h2, xMin, xMax, Δϕ)  
	    numerator = r + sqrt(r)
	    denominator = 1-r
	    return nSigma^2*( numerator / denominator )^2
	end
	
	function get_M_min(r::Real, nSigma::Real)
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
end

# ╔═╡ 7d948d45-f76b-465f-be29-338d926b82a9
h1 = Hist1D(dfs[1].phi, 0:Δϕ:180  ); # reference histogram

# ╔═╡ a1fbd981-f4c6-4899-9156-be5629f5130c
h2 = Hist1D(dfs[2].phi, 0:Δϕ:180  ); # compared histogram

# ╔═╡ 6c8d381d-994a-45b3-b59e-73d3c3f4e29d
matRatios = zeros(Int(180/Δϕ),Int(180/Δϕ)); # container for each r_i

# ╔═╡ f1f65c2b-480e-463b-aa3d-ba551a097843
for ϕmin in 0:Δϕ:180-Δϕ        # iterating over the range ϕmin, ϕmax
    for ϕmax in ϕmin+Δϕ:Δϕ:180
        r = get_r(h1, h2, ϕmin, ϕmax, Δϕ)

		matRatios[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = r
	end
end

# ╔═╡ af03769b-34d1-420a-8397-709b8f6218a3
begin # this is just to clean up the zeros so they don't show in the map
	replace!( matRatios, 0.0 => NaN  ) ;
	replace!( matRatios, Inf => NaN  ) ;
	nothing
end

# ╔═╡ ef2bb53e-c19c-490d-a116-c3cba8a84c25
c   = :jet

# ╔═╡ efa8d8e0-8207-416f-8735-c82f8d48f0da
hmRatios = heatmap( 
    (Δϕ/2:Δϕ:180-Δϕ/2), 
    (Δϕ/2:Δϕ:180-Δϕ/2), 
    matRatios,
    xlabel = "ϕmin",
    ylabel = "ϕmax",
    colorbar_title= "r",
    title = "r; K: -0.65, - 1.0; Δϕ = $Δϕ",
    right_margin = 5Plots.mm,
    left_margin = 5Plots.mm,
    size = (800,600),
    c   = c, 
)

# ╔═╡ 3be41328-baa6-40f5-b76f-1e5e8c31258c
begin
	minR = 1e6
	idx = (0,0)
	for i in CartesianIndices(matRatios)
		if( matRatios[i] != NaN && matRatios[i] < minR )
			minR = matRatios[i]
			idx =  i[1]*Δϕ-Δϕ, i[2]*Δϕ
		end
	end
	minR, idx = round(minR, digits = 2), idx
end

# ╔═╡ a5037ddd-09bd-4f59-b9b8-c23dba5de9b8
md"""
In the figure we can see the calculated ``r_i`` for each range of ``\phi``'s. The most ideal value is the lowest, ie. r = $minR for ``\phi \in`` $idx ``{}^{\circ}``. 

This however, does not tell the whole story yet. We must consider the uncertainty of the result as well as the uncertainty of the simulation. Furthermore, since the aim was to answer **how many events are required** to distinguish the two spectra, we still have some steps to take. 

We will therefore produce a few more figures and calculations. 
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

Lastly, since the behavior of r is very much non-linear (``\tilde M`` rapidly explodes for ``r`` close to 1.0) we consider the maximum (and minimum) uncertainties:

``r_{max} = \frac{M + \Delta M }{N - \Delta N}``
``r_{min} = \frac{M - \Delta M}{N - \Delta N}``.

Here, the worst case scenario (``r_{max}``) can in some cases exceed 1.0, (ie. ``r > 1.0``). Such regions will be removed from the analysis and will show in the maps as empty squares. 

Finally, let us look at all of the mentioned values. 
In the figure below we show maps for: 1. ``r``, 2.``\Delta r``, 3. ``\Delta r / r *100``, 4. ``r_{min}``, 5. ``r_{max}``, 6. ``\tilde M(r)``, 7. ``\tilde M(r_{min})``, and 8. ``\tilde M(r_{max})``. 
"""

# ╔═╡ 953e18ab-482a-4fc0-b993-8701a42b70ea
# initialize containers for the individual values
begin
	matDeltaRatios = zeros(Int(180/Δϕ),Int(180/Δϕ))
	matRelativeRatios = zeros(Int(180/Δϕ),Int(180/Δϕ))
	matRmin = zeros(Int(180/Δϕ),Int(180/Δϕ))
	matRmax = zeros(Int(180/Δϕ),Int(180/Δϕ))
	
	matMr = zeros(Int(180/Δϕ),Int(180/Δϕ))
	matMrmin = zeros(Int(180/Δϕ),Int(180/Δϕ))
	matMrmax = zeros(Int(180/Δϕ),Int(180/Δϕ))
	
	nothing
end

# ╔═╡ 2382ccbf-61c4-4c81-9c86-65b3d46ec082
for ϕmin in 0:Δϕ:180-Δϕ
    for ϕmax in ϕmin+Δϕ:Δϕ:180
        
        r = get_r(h1, h2, ϕmin, ϕmax, Δϕ) 
        Δr = get_delta_r( h1, h2, ϕmin, ϕmax, Δϕ )
        δr = Δr/r*100
        rmin = get_r_min(h1, h2, ϕmin, ϕmax, Δϕ) 
        rmax = get_r_max(h1, h2, ϕmin, ϕmax, Δϕ) <= 1.0 ? get_r_max(h1, h2, ϕmin, ϕmax, Δϕ) : NaN
    
        matDeltaRatios[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = Δr
        matRelativeRatios[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = δr
        matRmin[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = rmin
        matRmax[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = rmax
        
        matMr[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = rmax <= 1.0 ? get_M_min( r, 1.0 ) : NaN
        matMrmin[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = rmax <= 1.0 ? get_M_min( rmin, 1.0 ) : NaN
        matMrmax[Int(ϕmax/Δϕ), Int(ϕmin/Δϕ)+1] = get_M_min( rmax, 1.0 )
    end
end

# ╔═╡ 5d0e2fa5-7a39-4671-b64f-07c0ee650eba
# cleaning up
begin
	replace!( matDeltaRatios, 0.0 => NaN  ) 
	replace!( matDeltaRatios, Inf => NaN  ) 
	
	replace!( matRelativeRatios, 0.0 => NaN  ) 
	replace!( matRelativeRatios, Inf => NaN  ) 
	
	replace!( matRmin, 0.0 => NaN  ) 
	replace!( matRmin, Inf => NaN  )
	
	replace!( matRmax, 0.0 => NaN  ) 
	replace!( matRmax, Inf => NaN  )
	
	
	replace!( matMr, 0.0 => NaN  ) 
	replace!( matMr, Inf => NaN  )
	
	replace!( matMrmin, 0.0 => NaN  ) 
	replace!( matMrmin, Inf => NaN  )
	
	replace!( matMrmax, 0.0 => NaN  ) 
	replace!( matMrmax, Inf => NaN  );
	
	nothing
end

# ╔═╡ 6c7072f3-5016-4cd9-ba99-a98cccb41dd2
begin
	 
	hmRelativeRatios = heatmap( 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    matRelativeRatios,
	    xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "Δr/r*100",
	    title = "Δr/r*100; κ: -0.65 - 1.0; Δϕ = $Δϕ",
	    right_margin = 5Plots.mm,
	    size = (800,600),
	    c   = c, 
	)
	
	
	hmRmin = heatmap( 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    matRmin,
	    xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "rMin ",
	    title = "rMin; κ: -0.65 - 1.0; Δϕ = $Δϕ",
	    right_margin = 10Plots.mm,
	    size = (800,600),
	    c   = c, 
	)
	
	hmRmax = heatmap( 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    matRmax,
	    xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "rMax ",
	    title = "rMax; κ: -0.65 - 1.0; Δϕ = $Δϕ",
	    right_margin = 10Plots.mm,
	    size = (800,600),
	    c   = c, 
	)
	
	
	
	hmMr = heatmap( 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    matMr ,
	    xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "M(r) ",
	    title = "M(r); κ: -0.65 - 1.0; Δϕ = $Δϕ",
	    right_margin = 10Plots.mm,
	    size = (800,600),
	    c   = c, 
	    colorbar_scale = :log10,
	)
	
	hmMrmin = heatmap( 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    matMrmin ,
	    xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "M(rmin) ",
	    title = "M(rmin); κ: -0.65 - 1.0; Δϕ = $Δϕ",
	    right_margin = 10Plots.mm,
	    size = (800,600),
	    c   = c, 
	    colorbar_scale = :log10,
	)
	
	hmMrmax = heatmap( 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    (Δϕ/2:Δϕ:180-Δϕ/2), 
	    matMrmax,
	    xlabel = "ϕmin",
	    ylabel = "ϕmax",
	    colorbar_title= "M(rmax)",
	    title = "M(rmax); κ: -0.65 - 1.0; Δϕ = $Δϕ",
	    right_margin = 10Plots.mm,
	    size = (800,600),
	    c   = c, 
	    colorbar_scale = :log10,
	)

	nothing
end

# ╔═╡ eba85150-3ec9-48aa-a441-2c79204d8ab3
plot(
    hmRatios, 
    hmRelativeRatios,  
    hmRelativeRatios,
    hmRmin, 
    hmRmax, 
    hmMr,
    hmMrmin, 
    hmMrmax,
    left_margin= 10Plots.mm,
    right_margin= 10Plots.mm,
    size =(2200, 1600) 
)

# ╔═╡ d322dbc1-7dbc-436b-99b4-d90f424fefe5
begin
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
			
end

# ╔═╡ 08403fbf-edef-4212-b8bf-57547e8c8975
begin
	minRrelativeRatio = round(get_min_value(matRelativeRatios), digits = 2)
	idxRrelativeRatio = get_min_idx(matRelativeRatios)[1]*Δϕ-Δϕ, get_min_idx(matRelativeRatios)[2]*Δϕ

	minMr = round(get_min_value(matMr), digits = 2)
	idxMr = get_min_idx(matMr)[1]*Δϕ-Δϕ, get_min_idx(matMr)[2]*Δϕ
	
	minMmin = round(get_min_value(matMrmin), digits = 2)
	idxMmin = get_min_idx(matMrmin)[1]*Δϕ-Δϕ, get_min_idx(matMrmin)[2]*Δϕ

	minMmax = round(get_min_value(matMrmax), digits = 2)
	idxMmax = get_min_idx(matMrmax)[1]*Δϕ-Δϕ, get_min_idx(matMrmax)[2]*Δϕ

end

# ╔═╡ b0631780-5d84-4033-a920-cfde28008aa8
md"""
The resultant minimal values are:

 + ``r`` = $minR @ $idx``{}^{\circ}`` 
 + ``\delta r`` = $minRrelativeRatio @ $idxRrelativeRatio``{}^{\circ}``
 + ``M(r)`` = $minMr @ $idxMr``{}^{\circ}``
 + ``M(r_{min})`` = $minMmin @ $idxMmin``{}^{\circ}``
 + ``M(r_{max})`` = $minMmax @ $idxMmax``{}^{\circ}``
"""

# ╔═╡ ed4d68ee-131a-4a44-99a3-eac9ece030c8
md"""
This is as far as we got so far. Conclusions are not yet made. Further analysis will be done later. 
"""

# ╔═╡ Cell order:
# ╟─afeccb00-7abf-11ed-194a-558eaaa68040
# ╟─411eeb6b-6902-4cd9-8bda-8ebae1de81b1
# ╟─5565017c-534f-46a5-b86b-37753c42da24
# ╟─68624ab0-a258-4132-8e5a-54da853d856d
# ╟─66154af6-bb40-46d8-84c4-3e350acb8676
# ╟─a7e6341c-6f8b-42a6-b440-368fea258855
# ╟─6424e544-8aae-4e2e-82ce-f07c0c9e8083
# ╟─cad5c29a-f701-478d-902a-1a523137b5ee
# ╟─6c4cc9cc-3344-4ff0-bed3-df7ba1d00bc7
# ╟─739884af-76d1-46e8-8ee4-332c9baf40b8
# ╟─7edf09c5-bf7d-4668-a546-279ba3bc23d7
# ╟─9afbcc00-97dc-4497-84a4-72a2a86cd16b
# ╟─5113de94-8b41-4f86-99fd-4382142edc28
# ╟─f06ff5ff-6890-482a-b461-bcb2d3c23fb5
# ╟─26eb6233-2fb5-471a-bcf7-bf6e014fbc0c
# ╟─e2595477-a045-4330-9ef9-a246be1a3112
# ╟─95a6b54a-a692-4477-954f-849f015e6825
# ╟─85cb8fd5-70c9-423a-8ba0-94a8632d5291
# ╟─b7d5abe0-a806-4a6d-acf5-c0f3b9a75d33
# ╟─7d948d45-f76b-465f-be29-338d926b82a9
# ╟─a1fbd981-f4c6-4899-9156-be5629f5130c
# ╟─6c8d381d-994a-45b3-b59e-73d3c3f4e29d
# ╟─f1f65c2b-480e-463b-aa3d-ba551a097843
# ╟─af03769b-34d1-420a-8397-709b8f6218a3
# ╟─ef2bb53e-c19c-490d-a116-c3cba8a84c25
# ╟─efa8d8e0-8207-416f-8735-c82f8d48f0da
# ╟─3be41328-baa6-40f5-b76f-1e5e8c31258c
# ╟─a5037ddd-09bd-4f59-b9b8-c23dba5de9b8
# ╟─31f648a2-cdb3-4def-8f04-dbfc65787536
# ╟─953e18ab-482a-4fc0-b993-8701a42b70ea
# ╟─2382ccbf-61c4-4c81-9c86-65b3d46ec082
# ╟─5d0e2fa5-7a39-4671-b64f-07c0ee650eba
# ╟─6c7072f3-5016-4cd9-ba99-a98cccb41dd2
# ╟─eba85150-3ec9-48aa-a441-2c79204d8ab3
# ╟─d322dbc1-7dbc-436b-99b4-d90f424fefe5
# ╟─08403fbf-edef-4212-b8bf-57547e8c8975
# ╟─b0631780-5d84-4033-a920-cfde28008aa8
# ╟─ed4d68ee-131a-4a44-99a3-eac9ece030c8
