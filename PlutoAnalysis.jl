### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 411eeb6b-6902-4cd9-8bda-8ebae1de81b1
# loading necessary packages
begin
    using Revise
    using StatsPlots,
        UnROOT,
        StatsBase,
        Polynomials,
        ColorSchemes,
        Suppressor,
        HypothesisTests,
        LaTeXStrings
    using FHist, DataFramesMeta, Distributions, DataFrames, RecipesBase
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
        size = (800, 800),
        legend = :topleft,
        guidefontsize = 16,
        tickfontsize = 12,
        titlefontsize = 16,
        legendfontsize = 12,
        left_margin = 4Plots.mm,
        right_margin = 8Plots.mm,
        top_margin = 4Plots.mm,
        bottom_margin = 6Plots.mm,
        dpi = 200,
        :colorbar_titlefontsize => 20,
        widen = :false,
        :markerstrokewidth => 1,
        :markerstrokecolor => :black,
    )

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
    lw = 4,
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

    file1 = ROOTFile(fName1)
    file2 = ROOTFile(fName2)

    singleElectronEnergies1 =
        AM.fill_from_root_file(file1, "tree", "reconstructedEnergy1") .+
        AM.fill_from_root_file(file1, "tree", "reconstructedEnergy2") # vector of single-electron energies for reference spectrum    
    singleElectronEnergies2 =
        AM.fill_from_root_file(file2, "tree", "reconstructedEnergy1") .+
        AM.fill_from_root_file(file2, "tree", "reconstructedEnergy2") # vector of single-electron energies for compared spectrum    

    phi1 = AM.fill_from_root_file(file1, "tree", "phi") # vector phi angles for reference spectrum    
    phi2 = AM.fill_from_root_file(file2, "tree", "phi") # vector phi angles for compared spectrum   
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
        normed = :true,
        nbins = 0:Δϕ:180,
        label = ["Κ=-0.88 (SM)" "Κ=-0.6639 (refined)"],
        xlabel = "escape angle φ [°]",
        ylabel = "normalized count rate",
        lw = 4,
    )

    resPhi = scatter(
        midpoints(0:Δϕ:180),
        AM.get_residuals(phi1, phi2, 0:Δϕ:180, true),
        yerr = AM.get_residuals_errors(phi1, phi2, 0:Δϕ:180),
        mc = :black,
        label = "",
        ylabel = "refined/SM",
        xlabel = "escape angle φ [°]",
        ms = 1.5,
        ylims = (0.8, 1.2),
    )

    shEne = histogram(
        [singleElectronEnergies1, singleElectronEnergies2],
        normed = :true,
        nbins = 0:ΔE:3000,
        label = ["SM" "ξ31=0.37 (refined)"],
        xlabel = L"\textrm{single ~electron ~energy ~T_e ~[keV]}",
        ylabel = "normalized count rate",
        legend = :best,
        lw = 4,
    )

    resEne = scatter(
        midpoints(0:ΔE:3000),
        AM.get_residuals(singleElectronEnergies1, singleElectronEnergies2, 0:ΔE:3000, true),
        yerr = AM.get_residuals_errors(
            singleElectronEnergies1,
            singleElectronEnergies2,
            0:ΔE:3000,
        ),
        mc = :black,
        label = "",
        ylabel = "refined/SM",
        xlabel = L"\textrm{single ~electron ~energy ~T_e ~[keV]}",
        ylims = (0.8, 1.2),
        ms = 1.5,
    )

    plot(
        shPhi,
        shEne,
        resPhi,
        resEne,
        plot_title = "reconstructed distributions",
        size = (2000, 1000),
        thickness_scaling = 1.4,
        layout = @layout [a b; c{0.3h} d]
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
plot(
    0.8:0.001:0.999,
    AM.Mmin.(0.8:0.001:0.999, 1e5, 1e5),
    xlabel = "r",
    ylabel = "M(r)",
    label = "",
    lw = 4,
    size = (800, 400),
    yscale = :log10,
    yticks = [10^3, 10^4, 10^5, 10^6, 10^7, 10^8],
    xticks = 0.8:0.02:1.0,
    thickness_scaling = 1.2,
)

# ╔═╡ 6c8d381d-994a-45b3-b59e-73d3c3f4e29d
begin
    BBBPhi = AM.BBB(phi1, phi2, 0.0, 180.0, float(Δϕ), 1.0)
    BBBEne = AM.BBB(
        singleElectronEnergies1,
        singleElectronEnergies2,
        0.0,
        3000.0,
        float(ΔE),
        1.0,
    )
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
        colorbar_title = "r",
        right_margin = 5Plots.mm,
        left_margin = 5Plots.mm,
        size = (800, 600),
        title = "r-ratios for BBBPhi",
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
        colorbar_title = "r",
        right_margin = 5Plots.mm,
        left_margin = 5Plots.mm,
        size = (800, 600),
        title = "r-ratios for BBBEne",
    )

    plot(
        hmrBBBPhi,
        hmrBBBEne,
        aspect_ratio = 1,
        size = (1600, 750),
        thickness_scaling = 1.4,
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
        colorbar_title = "\nStats needed",
        right_margin = 5Plots.mm,
        left_margin = 5Plots.mm,
        size = (800, 600),
        title = "Stats needed for BBBPhi",
        colorbar_scale = :log10,
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
        colorbar_title = "\nStats needed",
        right_margin = 15Plots.mm,
        left_margin = 5Plots.mm,
        size = (800, 600),
        title = "Stats needed for BBBEne",
        colorbar_scale = :log10,
    )

    plot(hmsBBBPhi, hmsBBBEne, size = (1600, 700), thickness_scaling = 1.4)
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
    sampleSizes = vcat(collect(20_000:10_000:100_000), collect(150_000:50_000:800_000))
    xticks = (1:length(sampleSizes), sampleSizes)
    CL = 0.95
end

# ╔═╡ a593dd52-19e8-4502-9ff5-2211c6a4565a
begin
	Chi2Ene = AM.Chi2(
		singleElectronEnergies1, singleElectronEnergies2, 450, 2100, 150, sampleSizes,  CL)
	Chi2Phi = AM.Chi2(phi1, phi2, 0, 180, 15, sampleSizes,  CL)
end

# ╔═╡ 576f03a1-1d95-4a0a-a8df-d368af5d6d37
begin
	KSEne = @suppress AM.KS(singleElectronEnergies1, singleElectronEnergies2,  CL)
	KSPhi = @suppress AM.KS(phi1, phi2, CL)
end

# ╔═╡ b70e1b8d-935c-4625-ace4-347924a00570
begin
	pValsKSPhi = AM.get_pVals(KSPhi, sampleSizes)
	pValsKSEne = AM.get_pVals(KSEne, sampleSizes)

	pValsChi2Phi = Chi2Phi.pVals
	pValsChi2Ene = Chi2Ene.pVals
end

# ╔═╡ c08eec73-44df-42c2-845e-5f654093cddc
begin
	meansKSPhi = mean.(pValsKSPhi)
	meansKSEne = mean.(pValsKSPhi)
	meansChi2Phi = mean.(pValsChi2Phi)
	meansChi2Ene = mean.(pValsChi2Ene)
end


# ╔═╡ e7e5a693-efe1-4a2d-8aa1-4c24ee377168
md"""
**Disclaimer: For ``\chi^2`` there is an extra condition that must be filfilled taking into account the methodology. The binning must be chosen so that at least 5 events are in each bin. If this is not fulfilled, it will throw error.**
"""

# ╔═╡ dbbe2190-9af5-43c5-b9b5-6f5d5cab77c2
plot(
	[meansKSPhi, meansKSEne, meansChi2Phi, meansChi2Ene],
	label = ["KS: φ" "KS: Energy" "χ2: φ" "χ2: Energy"],
	markersize = 8, 
	markershape = [:c :diamond :d :star5],
	xticks = xticks,
	xlabel = "sample size", 
	ylabel ="p-value",
	title= "mean p-values for combination of tests and distributions",
	xrotation = 55
)

# ╔═╡ 852b33c3-99fb-4d3f-b22c-737a301a9cd8
md"""
Above is a plot of mean p-values for each combination. We can see that the `energy` distribution needs more statistics to reject ``H_0``, in fact we have not even reach enough to obtain desired statistics. We look a bit closer at the data by looking at a boxplot, which shows the median and 2nd, 3rd quantile and outliers. 
"""

# ╔═╡ bb61aa57-0649-48fc-99e6-859c05f3e7a8
begin
	bxPhiKS = boxplot(pValsKSPhi, label = "", xticks=(1:length(sampleSizes), sampleSizes), xrotation = 45, c=1, fa = 0.5, xlabel = "sample size", ylabel ="p-value", title= "φ: KS test for various sample sizes", size =(1200,600), bottom_margin = 12Plots.mm, thickness_scaling = 1.5)

	
	lens!([length(sampleSizes)-12, length(sampleSizes)], [0.0, 0.05], inset = (1, bbox(0.57, 0.2, 0.4, 0.3)), lc=:black, lw =:3, xrotation=45, framestyle =:box)

	xaxis!(xticks=(length(sampleSizes)-12:length(sampleSizes), sampleSizes[end-12:end]), subplot =2)
end

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

# ╔═╡ 5e04ec18-cbcd-4487-858a-7b9dfd942e51
begin 
	bp1 = boxplot(
		pValsKSPhi, 
		label="", 
		xticks=:none, 
		xrotation=65, 
		c=1, 
		fa=0.5, 
		ylabel="p-value",
		top_margin=0Plots.px,
		bottom_margin=0Plots.px,
		left_margin=0Plots.px,
		right_margin=0Plots.px,
		ylims=(0,1),
		
	)
	annotate!([(18,0.8,("φ: KS", 20))])

	bp2 = boxplot(
		pValsKSEne, 
		label="", 
		xticks=:none, #(1:length(sampleSizes), sampleSizes), 
		xrotation=65, 
		c=2, 
		fa=0.5,  
		yticks=:none,
		top_margin=0Plots.px,
		bottom_margin=0Plots.px,
		left_margin=0Plots.px,
		right_margin=0Plots.px,
		ylims=(0,1)
	)
	annotate!([(18,0.8,("E: KS", 20))])

	bp3 = boxplot(
		filter!.(x -> x .>=0, pValsChi2Phi), 
		label="", 
		xticks=(1:length(sampleSizes), sampleSizes), 
		xrotation=65, 
		c=3, 
		fa=0.5, 
		xlabel="sample size", 
		ylabel ="p-value", 
		top_margin=0Plots.px,
		left_margin=0Plots.px,
		right_margin=0Plots.px,
		ylims=(0,1)
	)
	annotate!([(18,0.8,(L"\mathrm{φ:~ \chi^2}", 20))])
	

	bp4 = boxplot(
		filter!.(x -> x .>=0, pValsChi2Ene), 
		label="", 
		xticks=(1:length(sampleSizes), sampleSizes), 
		xrotation=65, 
		c=4, 
		fa=0.5, 
		xlabel="sample size", 
		top_margin=0Plots.px,
		left_margin=0Plots.px,
		right_margin=0Plots.px, 
		yticks=:false,
		ylims=(0,1)
	)
	annotate!([(18,0.8,(L"\mathrm{E:~ \chi^2}", 20))])
	

	bxpAll = plot(
		bp1, bp2, bp3, bp4, 
		size=(1200, 1100),
		thickness_scaling=1.3, 
		grid=:false, 
		minorgrid=:false
	)

end


# ╔═╡ 75f2bc28-529e-47e1-95cf-e8a6b2f864ec
md"""
However, since it can be seen that the p-values vary greatly, we can perform an additional analysis by looking at the efficiency of rejecting H0 for various sample sizes. That is, for each subset of size ``N`` we compute efficiency as ``\varepsilon = \frac{N_{reject}}{N}``. 
"""

# ╔═╡ dc755feb-58eb-4be2-a6f8-5e2b51a06411
begin
	alpha = 1-0.95
	
	effKSPhi95 = 
		count.(p -> p<=alpha, pValsKSPhi)./length.(pValsKSPhi).*100 .|> round .|>Int
	effKSEne95 = 
		count.(p -> p<=alpha, pValsKSEne)./length.(pValsKSEne).*100 .|> round .|>Int 
	
	effChiPhi95 = 
		count.(p -> p<=alpha, pValsChi2Phi)./length.(pValsChi2Phi).*100 .|> round .|>Int
	effChiEne95 = 
		count.(p -> p<=alpha, pValsChi2Ene)./length.(pValsChi2Ene).*100 .|> round .|>Int

end

# ╔═╡ 3651c0a2-acd0-4193-9cc4-d3f7d1cd55b3
begin
	grpbKS = groupedbar(
		names = repeat(1:length(sampleSizes), outer = 2),
		hcat(effKSPhi95, effKSEne95),
		groups = repeat( ["KS: φ", "KS: E"], inner = length(sampleSizes) ),
		c = [2 1],
		yerr = vcat(errKSPhi, errKSEne),
		ylabel = "efficiency", 
	    xlabel = "sample size",
		title  = L"\textrm{efficiency ~of ~passing ~KS}~\textrm{test ~for ~CL ~=} %$(CL)",
		xticks = (1:length(sampleSizes), sampleSizes),
		xrotation=45,
	)
	
	grpbChi = groupedbar(
		names = repeat(1:length(sampleSizes), outer = 2),
		hcat(effChiPhi95, effChiEne95),
		groups = repeat( ["Chi: φ", "Chi: E"], inner = length(sampleSizes) ),
		c = [4 3],
		yerr = vcat(errChiPhi, errChiEne),
		ylabel = "efficiency", 
	    xlabel = "sample size",
		title  = L"\textrm{efficiency ~of ~passing }~\mathrm{\chi^2} ~\textrm{test ~for ~CL ~=} %$(CL)",
		xticks = (1:length(sampleSizes), sampleSizes),
		xrotation=45,
	)
	
	plot(grpbKS, grpbChi, layout = grid(2,1))
end

# ╔═╡ 5163e320-0711-42e8-89c6-90994f612c56
md"""
Finally we can create a table of *S* for the various defined methods and sample sizes.
"""

# ╔═╡ 8a705628-67f8-45c1-be41-a0531e5378b5
for (cl, nsig) in zip([0.68, 0.90, 0.95], [1, 1.645, 1.96])
	minKSPhi = AM.get_best_sample_size(KSPhi, sampleSizes, cl, 100,3)
	minKSEne = AM.get_best_sample_size(KSEne, sampleSizes, cl, 100,3)
	minChi2Phi = AM.get_best_sample_size(Chi2Phi, sampleSizes, cl, 100,3)
	minChi2Ene = AM.get_best_sample_size(Chi2Ene, sampleSizes, cl, 100,3)
	minBBBPhi = AM.BBB(phi1, phi2, 0, 180, Δϕ, nsig).minEvents
	minBBBEne = AM.BBB(singleElectronEnergies1, singleElectronEnergies2, 0, 3000, ΔE, nsig).minEvents
end

# ╔═╡ 15149b75-4ce7-472c-adfa-8d4ad153ba95
md"""
CL\method | KSphi |KSEne | Χ²Phi | Χ²Ene | BBBPhi | BBB Ene
----------|-------|------|-------|-------|--------|--------
68%       |70k    |450k  |60k    |350k   |11k     |60k
90%       |150k   |750k  |60k    |350k   |30k     |160k
95%       |200k   |750k  |200k   |500k   |45k     |230k

"""

# ╔═╡ 9a2dbe13-66c9-4fa5-8404-2872601bafef

function get_delta_r_map(
    vec1::Vector{<:Real},
    vec2::Vector{<:Real},
    xMin::Real,
    xMax::Real,
    stepSize::Real,
)
    h1 = Hist1D(vec1, xMin:stepSize:xMax)
    h2 = Hist1D(vec2, xMin:stepSize:xMax)

    matRatios = zeros(Int(xMax / stepSize), Int(xMax / stepSize))

    for minROI = xMin:stepSize:xMax-stepSize        # iterating over the range xMin, xMax
        for maxROI = minROI+stepSize:stepSize:xMax
            rmin = AM.get_delta_r(h1, h2, minROI, maxROI, stepSize)

            matRatios[Int(maxROI / stepSize), Int(minROI / stepSize)+1] = rmin
        end
    end

    return matRatios
end



# ╔═╡ b675e1b6-320c-4db5-b527-068e6e571c84
get_delta_r_map(bbb::AM.BBB) =
    get_delta_r_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)

# ╔═╡ 45b56c40-2c2b-4725-9cf8-85958947f25b
get_r_min_map(bbb::AM.BBB) =
    AM.get_r_map(bbb.vector1, bbb.vector2, bbb.xMin, bbb.xMax, bbb.stepSize)

# ╔═╡ 9ee70dbc-e29a-4bfc-a5f8-84afe88eb00b
figDir = "/media/shoram/Extra SSD/CernBox/Work/Presentations/2023/2023_05_11_AnalysisMeeting/figs"

# ╔═╡ 02e4fab7-c042-4f80-a616-2f318f177302
begin
    a = get_r_min_map(BBBPhi)
    replace!(a, 0.0 => NaN)
    title = "r_min"
    afigtitle = "r_min_phi.png"

    b = get_r_min_map(BBBEne)
    replace!(b, 0.0 => NaN)
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
        size = (800, 600),
        title = title,
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
        size = (800, 600),
        title = title,
        colorbar_scale = :log10,
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

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "ff9bffc15194090e4a842a319e3e9e23ff75cc69"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8bc0aaec0ca548eb6cf5f0d7d16351650c1ee956"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.2"
weakdeps = ["ChainRulesCore"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

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

[[deps.BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "fc54d5837033a170f3bad307f993e156eefc345f"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.2.7"

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
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "42fe66dbc8f1d09a44aa87f18d26926d06a35f84"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.3"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "d730914ef30a06732bdd9f763f6cc32e92ffbff1"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.3.1"

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
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

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
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "4e88377ae7ebeaf29a047aa1ee40826e0b708a5d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.7.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "db40d3aff76ea6a3619fdd15a8c78299221a2394"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.97"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"

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

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

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
git-tree-sha1 = "b4fbdd20c889804969571cc589900803edda16b7"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.7.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FHist]]
deps = ["BayesHistogram", "LinearAlgebra", "MakieCore", "Measurements", "RecipesBase", "Requires", "Statistics", "StatsBase", "UnicodePlots"]
git-tree-sha1 = "445c9b5da795962ccfa7747daae5b05a4902ff28"
uuid = "68837c9b-b678-4cd5-9925-8a54edc8f695"
version = "0.9.3"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "0b3b52afd0f87b0a3f5ada0466352d125c9db458"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.2.1"

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
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8b8a2fd4536ece6e554168c21860b6820a8a83db"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "19fad9cd9ae44847fe842558a744748084a722d1"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.7+0"

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
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "2613d054b0e18a3dea99ca1594e9a3960e025da4"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.7"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "0ec02c648befc2f94156eaef13b0f38106212f3f"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.17"

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
git-tree-sha1 = "fee0691e3336a71503dada09ed61bed786b0f59f"
uuid = "09f84164-cd44-5f33-b23f-e6b0d136a0d5"
version = "0.10.13"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "d75853a0bdbfb1ac815478bacd89cd27b550ace6"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.3"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0cb9352ef2e01574eeebdb102948a58740dcaf83"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2023.1.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IterTools]]
git-tree-sha1 = "4ced6667f9974fc5c5943fa5e2ef1ca43ea9e450"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.8.0"

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
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "6a125e6a4cb391e0b9adbd1afa9e771c2179f8ef"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.23"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "90442c50e202a5cdf21a7899c66b240fdef14035"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.7"

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

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LRUCache]]
git-tree-sha1 = "48c10e3cc27e30de82463c27bef0b8bdbd1dc634"
uuid = "8ac3fa9e-de4c-5943-b1dc-09c6b5f20637"
version = "1.4.1"

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
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

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
git-tree-sha1 = "1258d920c598f9d5a5378284aed33a96e991c974"
uuid = "9255714d-24a7-4b30-8ea3-d46a97f7e13b"
version = "0.4.2"

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
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "c3ce8e7420b3a6e071e0fe4745f5d4300e37b13f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.24"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.LorentzVectors]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "da81cd93228bfa12d8d9dd8aaab61f48e5f74441"
uuid = "3f54b04b-17fc-5cd4-9758-90c048d965e3"
version = "0.4.3"

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
git-tree-sha1 = "154d7aaa82d24db6d8f7e4ffcfe596f40bff214b"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2023.1.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9926529455a331ed73c19ff06d16906737a876ed"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.3"

[[deps.MarchingCubes]]
deps = ["PrecompileTools", "StaticArrays"]
git-tree-sha1 = "c8e29e2bacb98c9b6f10445227a8b0402f2f173a"
uuid = "299715c1-40a9-479a-aaf9-4a633d36f717"
version = "0.1.8"

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
version = "2.28.2+0"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "51d946d38d62709d6a2d37ea9bcc30c80c686801"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.9.0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Memoization]]
deps = ["MacroTools"]
git-tree-sha1 = "073f080e733bc6697411901224ed4fd15fefaffa"
uuid = "6fafb56a-5788-4b4e-91ca-c0cea6611c73"
version = "0.2.1"

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
version = "2022.10.11"

[[deps.MultivariateStats]]
deps = ["Arpack", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "68bf5103e002c44adfd71fea6bd770b3f0586843"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.2"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

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
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa4b74f80b01c6bc2b89992b861b5f210e665b5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.21+0"

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
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "4b2e829ee66d4218e0cef22c0a64ee37cf258c29"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "75ca67b2c6512ad2d0c767a7cfc55e75075f8bbc"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.16"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "b478a748be27bd2f2c73a7690da219d0844db305"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.51"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase"]
git-tree-sha1 = "3aa2bb4982e575acd7583f01531f241af077b163"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.2.13"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "213579618ec1f42dea7dd637a42785a608b1ea9c"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.4"

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
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

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
git-tree-sha1 = "1e597b93700fa4045d7189afa7c004e0584ea548"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.3"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "de432823e8aab4dd1a985be4be768f95acf152d4"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.17"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "04bdff0b09c65ff3e06a05e3eb7b120223da3d39"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.0"

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
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "7beb031cf8145577fbccacd94b8a8f4ce78428d3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "832afbae2a45b4ae7e831f86965469a24d1d8a83"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.26"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "14ef622cf28b05e38f8af1de57bc9142b03fbfe3"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.5"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.Suppressor]]
deps = ["Logging"]
git-tree-sha1 = "37d1976ca8368f6adbe1d65a4deeeda6ee7faa31"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.4"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

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
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.UnROOT]]
deps = ["AbstractTrees", "ArraysOfArrays", "BitIntegers", "CodecLz4", "CodecXz", "CodecZstd", "HTTP", "IterTools", "LRUCache", "LibDeflate", "LorentzVectors", "Memoization", "Mixers", "Mmap", "Parameters", "PrettyTables", "SentinelArrays", "StaticArrays", "StructArrays", "Tables", "xrootdgo_jll"]
git-tree-sha1 = "f592e6a781f27e848210d718ed640c00e613c1ad"
uuid = "3cd96dde-e98d-4713-81e9-a4a1b0235ce9"
version = "0.8.23"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.UnicodePlots]]
deps = ["ColorSchemes", "ColorTypes", "Contour", "Crayons", "Dates", "LinearAlgebra", "MarchingCubes", "NaNMath", "PrecompileTools", "Printf", "Requires", "SparseArrays", "StaticArrays", "StatsBase"]
git-tree-sha1 = "b96de03092fe4b18ac7e4786bee55578d4b75ae8"
uuid = "b8865327-cd53-5732-bb35-84acbb429228"
version = "3.6.0"

    [deps.UnicodePlots.extensions]
    FreeTypeExt = ["FileIO", "FreeType"]
    ImageInTerminalExt = "ImageInTerminal"
    IntervalSetsExt = "IntervalSets"
    TermExt = "Term"
    UnitfulExt = "Unitful"

    [deps.UnicodePlots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    FreeType = "b38be410-82b0-50bf-ab77-7b57e271db43"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Term = "22787eb5-b846-44ae-b979-8e399b8463ab"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "ba4aa36b2d5c98d6ed1f149da916b3ba46527b2b"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.14.0"

    [deps.Unitful.extensions]
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8abe223c2549ea70be752b20a53aa236a7868eb0"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.3+0"

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
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

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
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libdeflate_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "42f15e41ea8cec5f460b4546170912f410c694f5"
uuid = "46979653-d7f6-5232-b59e-dd310c4598de"
version = "1.18.0+0"

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
# ╠═a593dd52-19e8-4502-9ff5-2211c6a4565a
# ╠═576f03a1-1d95-4a0a-a8df-d368af5d6d37
# ╠═b70e1b8d-935c-4625-ace4-347924a00570
# ╠═c08eec73-44df-42c2-845e-5f654093cddc
# ╟─e7e5a693-efe1-4a2d-8aa1-4c24ee377168
# ╠═dbbe2190-9af5-43c5-b9b5-6f5d5cab77c2
# ╟─852b33c3-99fb-4d3f-b22c-737a301a9cd8
# ╠═bb61aa57-0649-48fc-99e6-859c05f3e7a8
# ╟─51c2801e-7051-490a-bd68-1dceb9d7f7f9
# ╟─05700042-ca64-40d7-8ef2-76a75897fe11
# ╠═5e04ec18-cbcd-4487-858a-7b9dfd942e51
# ╟─75f2bc28-529e-47e1-95cf-e8a6b2f864ec
# ╠═dc755feb-58eb-4be2-a6f8-5e2b51a06411
# ╠═3651c0a2-acd0-4193-9cc4-d3f7d1cd55b3
# ╟─5163e320-0711-42e8-89c6-90994f612c56
# ╠═8a705628-67f8-45c1-be41-a0531e5378b5
# ╠═15149b75-4ce7-472c-adfa-8d4ad153ba95
# ╠═9a2dbe13-66c9-4fa5-8404-2872601bafef
# ╠═b675e1b6-320c-4db5-b527-068e6e571c84
# ╠═45b56c40-2c2b-4725-9cf8-85958947f25b
# ╠═02e4fab7-c042-4f80-a616-2f318f177302
# ╠═9ee70dbc-e29a-4bfc-a5f8-84afe88eb00b
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
