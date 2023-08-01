cd(@__DIR__) #where repository is at
using Pkg;
Pkg.activate(Base.current_project());


## Packages
using Revise
using StatsPlots, UnROOT, StatsBase, BenchmarkTools
using FHist, DataFramesMeta, Distributions, RecipesBase, Suppressor, HypothesisTests


## Including the Analysis Module 
push!(LOAD_PATH, pwd()*"/src/")
using AnalysisModule
Revise.track(AnalysisModule)


## Load data 
fName1 = "../Job17/Data/2vbb_Se82_Falaise_EneThetaPhi_1e8E.root"                           # root file for reference spectra
fName2 = "../Job17/Data/G0-G4_xi31_037_kappa_m06639_1e8E.root"   # root file for compared spectra

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

## Analysis 
CL = 0.50
sampleSizes = vcat(collect(20_000:10_000:100_000), collect(150_000:50_000:800_000))

Chi2Phi = Chi2(phi1, phi2, CL, 0, 180, 90, sampleSizes = sampleSizes)

KSPhi = KS(phi1, phi2, CL,sampleSizes = sampleSizes)

@benchmark KS(phi1, phi2, CL,sampleSizes = sampleSizes)