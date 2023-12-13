cd(@__DIR__) #where repository is at
using Pkg;
Pkg.activate(Base.current_project());


## Packages
using Revise
using StatsPlots, UnROOT, StatsBase, BenchmarkTools
using FHist, DataFramesMeta, Distributions, RecipesBase, Suppressor, HypothesisTests


## Including the Analysis Module 
push!(LOAD_PATH, pwd() * "/src/")
using AnalysisModule
Revise.track(AnalysisModule)


## Load data 
fName1 = joinpath("../Job21/Data_wo_Bfield/Falaise_2nubb_1e6E.root")                           # root file for reference spectra
fName2 = joinpath("../Job21/Data_wo_Bfield/Xi31_060_1e6E.root")   # root file for compared spectra

file1 = ROOTFile(fName1)
file2 = ROOTFile(fName2)

singleElectronEnergies1 =
    vcat(
        fill_from_root_file(file1, "tree", "reconstructedEnergy1"),
        fill_from_root_file(file1, "tree", "reconstructedEnergy2")
    ) # vector of single-electron energies for reference spectrum    
singleElectronEnergies2 =
    vcat(
        fill_from_root_file(file2, "tree", "reconstructedEnergy1"),
        fill_from_root_file(file2, "tree", "reconstructedEnergy2")
    ) # vector of single-electron energies for compared spectrum    

phi1 = fill_from_root_file(file1, "tree", "phi") # vector phi angles for reference spectrum    
phi2 = fill_from_root_file(file2, "tree", "phi") # vector phi angles for compared spectrum    

## Analysis 
CL = 0.9
sampleSizes = vcat(collect(10_000:10_000:100_000), collect(150_000:50_000:500_000))

sh1 = stephist(phi1, nbins = 0:5:180, normed = :true, subplot = 1)
stephist!(phi2, nbins = 0:5:180, normed = :true)

sh2 = stephist(singleElectronEnergies1, nbins = 50:50:1900, normed = :true)
stephist!(singleElectronEnergies2, nbins = 50:50:1900, normed = :true)

# plot(sh1, sh2)
Chi2Phi = Chi2(phi1, phi2, CL, 0, 180, 15; sampleSizes=sampleSizes)
# KSPhi = KS(phi1, phi2, CL; sampleSizes=sampleSizes)

Chi2Ene = Chi2(singleElectronEnergies1, singleElectronEnergies2, CL, 50, 1200, 50; sampleSizes=sampleSizes)
KSEne = KS(singleElectronEnergies1, singleElectronEnergies2, CL; sampleSizes=sampleSizes)


histogram(fill_from_root_file(file1, "tree", "reconstructedEnergy1"), bins = 0:150:3000, normed = :true, label = "first")
histogram!(fill_from_root_file(file1, "tree", "reconstructedEnergy2"), bins = 0:150:3000, normed = :true, label = "second")


sh2 = stephist(singleElectronEnergies1[1:20_000], nbins = 50:50:1700, normed = :true)
stephist!(singleElectronEnergies2[1:20_000], nbins = 50:50:1700, normed = :true)

