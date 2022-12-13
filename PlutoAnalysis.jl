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


# ╔═╡ 953e18ab-482a-4fc0-b993-8701a42b70ea


# ╔═╡ 2382ccbf-61c4-4c81-9c86-65b3d46ec082


# ╔═╡ 5d0e2fa5-7a39-4671-b64f-07c0ee650eba


# ╔═╡ 6c7072f3-5016-4cd9-ba99-a98cccb41dd2


# ╔═╡ eba85150-3ec9-48aa-a441-2c79204d8ab3


# ╔═╡ ed4d68ee-131a-4a44-99a3-eac9ece030c8


# ╔═╡ 3ac14307-ec19-445a-91b7-1e798dd67b28


# ╔═╡ f02a21e9-8676-4b08-9a0f-01775d5e45ff


# ╔═╡ 71bd6ca6-7b5a-41ca-ab6d-c51f3d151c03


# ╔═╡ 4e3f5697-3461-4975-bed9-ffe4dfa2cc94


# ╔═╡ 9f68cd3e-0493-4b63-8369-2d2c2aa2e2c7


# ╔═╡ 411cb889-a4eb-4a87-a710-28cd3b89b36e


# ╔═╡ 7c637ea8-20fc-47f1-9026-5074a7039b93


# ╔═╡ e5d8955e-d482-4d14-b31b-c8d69bf12fed


# ╔═╡ 986cf167-96d3-437c-b59b-7a0fb0c0de77


# ╔═╡ 6fa583e3-1634-407f-8a94-93504c551bf2


# ╔═╡ 68abe250-d0e9-4c65-8a32-2c95b20397a4


# ╔═╡ 843af440-998e-4cbf-99da-10898237f8d3


# ╔═╡ 3a8bb8d9-4eff-4d66-9f90-65aa08e5145e


# ╔═╡ 99443bf5-a713-4bec-9dd6-889a9aa83972


# ╔═╡ 5481cb52-858c-476d-9cd6-69c0ac8ce352


# ╔═╡ a3652b27-bb87-4912-a87f-a752084a25fb


# ╔═╡ d40494d7-5a72-4f0f-806f-37cc63bc89d6


# ╔═╡ Cell order:
# ╟─afeccb00-7abf-11ed-194a-558eaaa68040
# ╠═411eeb6b-6902-4cd9-8bda-8ebae1de81b1
# ╠═5565017c-534f-46a5-b86b-37753c42da24
# ╠═68624ab0-a258-4132-8e5a-54da853d856d
# ╠═66154af6-bb40-46d8-84c4-3e350acb8676
# ╠═a7e6341c-6f8b-42a6-b440-368fea258855
# ╠═6424e544-8aae-4e2e-82ce-f07c0c9e8083
# ╠═cad5c29a-f701-478d-902a-1a523137b5ee
# ╠═6c4cc9cc-3344-4ff0-bed3-df7ba1d00bc7
# ╠═739884af-76d1-46e8-8ee4-332c9baf40b8
# ╠═7edf09c5-bf7d-4668-a546-279ba3bc23d7
# ╟─9afbcc00-97dc-4497-84a4-72a2a86cd16b
# ╠═5113de94-8b41-4f86-99fd-4382142edc28
# ╠═f06ff5ff-6890-482a-b461-bcb2d3c23fb5
# ╠═26eb6233-2fb5-471a-bcf7-bf6e014fbc0c
# ╠═e2595477-a045-4330-9ef9-a246be1a3112
# ╟─95a6b54a-a692-4477-954f-849f015e6825
# ╠═85cb8fd5-70c9-423a-8ba0-94a8632d5291
# ╠═953e18ab-482a-4fc0-b993-8701a42b70ea
# ╠═2382ccbf-61c4-4c81-9c86-65b3d46ec082
# ╠═5d0e2fa5-7a39-4671-b64f-07c0ee650eba
# ╠═6c7072f3-5016-4cd9-ba99-a98cccb41dd2
# ╠═eba85150-3ec9-48aa-a441-2c79204d8ab3
# ╠═ed4d68ee-131a-4a44-99a3-eac9ece030c8
# ╠═3ac14307-ec19-445a-91b7-1e798dd67b28
# ╠═f02a21e9-8676-4b08-9a0f-01775d5e45ff
# ╠═71bd6ca6-7b5a-41ca-ab6d-c51f3d151c03
# ╠═4e3f5697-3461-4975-bed9-ffe4dfa2cc94
# ╠═9f68cd3e-0493-4b63-8369-2d2c2aa2e2c7
# ╠═411cb889-a4eb-4a87-a710-28cd3b89b36e
# ╠═7c637ea8-20fc-47f1-9026-5074a7039b93
# ╠═e5d8955e-d482-4d14-b31b-c8d69bf12fed
# ╠═986cf167-96d3-437c-b59b-7a0fb0c0de77
# ╠═6fa583e3-1634-407f-8a94-93504c551bf2
# ╠═68abe250-d0e9-4c65-8a32-2c95b20397a4
# ╠═843af440-998e-4cbf-99da-10898237f8d3
# ╠═3a8bb8d9-4eff-4d66-9f90-65aa08e5145e
# ╠═99443bf5-a713-4bec-9dd6-889a9aa83972
# ╠═5481cb52-858c-476d-9cd6-69c0ac8ce352
# ╠═a3652b27-bb87-4912-a87f-a752084a25fb
# ╠═d40494d7-5a72-4f0f-806f-37cc63bc89d6
