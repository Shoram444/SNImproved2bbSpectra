using Pkg;
Pkg.activate(Base.current_project());


## Packages
using Revise
using StatsPlots, UnROOT, StatsBase, BenchmarkTools
using FHist, DataFramesMeta, Distributions, RecipesBase, Suppressor, HypothesisTests, LaTeXStrings


## Including the Analysis Module 
push!(LOAD_PATH, pwd() * "/src/")
using AnalysisModule
Revise.track(AnalysisModule)

fName1 = joinpath("../Job21/Data_wo_Bfield/Falaise_2nubb_1e6E.root")                           # root file for reference spectra
fName2 = joinpath("../Job21/Data_wo_Bfield/Xi31_060_1e6E.root")   # root file for compared spectra

file1 = ROOTFile(fName1)
file2 = ROOTFile(fName2)

singleElectronEnergies1 =
    fill_from_root_file(file1, "tree", "reconstructedEnergy1") .+
    fill_from_root_file(file1, "tree", "reconstructedEnergy2") # vector of single-electron energies for reference spectrum    
singleElectronEnergies2 =
    fill_from_root_file(file2, "tree", "reconstructedEnergy1") .+
    fill_from_root_file(file2, "tree", "reconstructedEnergy2") # vector of single-electron energies for compared spectrum    

phi1 = fill_from_root_file(file1, "tree", "phi") # vector phi angles for reference spectrum    
phi2 = fill_from_root_file(file2, "tree", "phi") # vector phi angles for compared spectrum  


stephist(
    singleElectronEnergies1, 
    nbins=0:50:3000, 
    size=(800, 500), 
    c=:black, lw=2, 
    label=L"2\nu\beta\beta ~standard", 
    xlabel=L"E_{sum} ~[\textrm{keV}]", 
    ylabel="normalized count rate", 
    normed=:true, 
    ylims=(0, 0.001),
    legend = :best,
    thickness_scaling = 1.5
)

stephist!(
    singleElectronEnergies2, 
    nbins=0:50:3000, 
    c=2, lw=2, 
    label=L"2\nu\beta\beta ~exotic", 
    normed = :true
)

savefig("Figs/Energy_comparison_standard_exotic.png")


stephist(
    phi1, 
    nbins=0:1:180, 
    size=(800, 500), 
    c=:black, lw=2, 
    label=L"2\nu\beta\beta ~standard", 
    xlabel=L"\theta ~[\degree]", 
    ylabel="normalized count rate", 
    normed=:true, 
    ylims=(0, 0.01),
    legend = :best,
    thickness_scaling = 1.5
)

stephist!(
    phi2, 
    nbins=0:1:180, 
    c=2, lw=2, 
    label=L"2\nu\beta\beta ~exotic", 
    normed = :true
)

savefig("Figs/Angular_comparison_standard_exotic.png")
