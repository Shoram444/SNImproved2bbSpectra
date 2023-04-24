# cd("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra")
using Pkg; Pkg.activate("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra")

using Revise 
using StatsPlots, UnROOT, StatsBase   
using FHist, DataFramesMeta, Distributions, RecipesBase, Suppressor

# include("/home/shoram/Work/PhD_Thesis/SNAngularCorrelation/AngularCorrelations/MiscFuncs.jl") # module with analysis functions 
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/AnalysisModule.jl")                # module with some extra miscelnacious functions
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/BBBModule.jl")                
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/Chi2Module.jl")              
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/KSModule.jl")         
include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/src/AnalysisModule.jl") 


# using .MiscFuncs, .AnalysisModule, .BBBModule, .Chi2Module, .KSModule
using .AnalysisModule

fName1 = "/home/shoram/Work/PhD_Thesis/Job17/Data/2vbb_Se82_Falaise_EneThetaPhi_1e8E.root"                           # root file for reference spectra
fName2 = "/home/shoram/Work/PhD_Thesis/Job17/Data/G0-G4_xi31_037_kappa_m06639_1e8E.root"   # root file for compared spectra

file1 = ROOTFile( fName1 )
file2 = ROOTFile( fName2 )

singleElectronEnergies1 = fill_from_root_file(file1, "tree","reconstructedEnergy1") .+ fill_from_root_file(file1, "tree", "reconstructedEnergy2") # vector of single-electron energies for reference spectrum    
singleElectronEnergies2 = fill_from_root_file(file2, "tree","reconstructedEnergy1") .+ fill_from_root_file(file2, "tree", "reconstructedEnergy2") # vector of single-electron energies for compared spectrum    

phi1 = fill_from_root_file(file1, "tree","phi") # vector phi angles for reference spectrum    
phi2 = fill_from_root_file(file2, "tree","phi") # vector phi angles for compared spectrum    

nSigma = 1
CL = 0.68
Δϕ = 15 # setting bin width
ΔE = 150

sampleSizes = vcat(collect(10_000:2000:18_000), collect(20_000:10_000:100_000), collect(100_000:50_000:400_000) )

BBBPhi = BBB(phi1, phi2, 0, 180, Δϕ, nSigma)
Chi2Phi = Chi2(phi1, phi2, 0, 180, Δϕ, sampleSizes,  CL)
KSPhi = @suppress KS(phi1, phi2, 0, 180, Δϕ, sampleSizes,  CL)

BBBEne = BBB(singleElectronEnergies1, singleElectronEnergies2, 0, 3000, ΔE, nSigma)
Chi2Ene = Chi2(singleElectronEnergies1, singleElectronEnergies2, 450, 2100, ΔE, sampleSizes,  CL)
KSEne = @suppress KS(singleElectronEnergies1, singleElectronEnergies2, 0, 3000, ΔE, sampleSizes,  CL)

@show BBBPhi.minEvents
@show Chi2Phi.minEvents
@show KSPhi.minEvents

@show BBBEne.minEvents
@show Chi2Ene.minEvents
@show KSEne.minEvents

