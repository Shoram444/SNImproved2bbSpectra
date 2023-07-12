# cd("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra")
using Pkg;
Pkg.activate("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra");

using Revise
using StatsPlots, UnROOT, StatsBase
using FHist, DataFramesMeta, Distributions, RecipesBase, Suppressor, HypothesisTests

# include("/home/shoram/Work/PhD_Thesis/SNAngularCorrelation/AngularCorrelations/MiscFuncs.jl") # module with analysis functions 
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/AnalysisModule.jl")                # module with some extra miscelnacious functions
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/BBBModule.jl")                
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/Chi2Module.jl")              
# include("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/KSModule.jl")         
includet("/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/src/AnalysisModule.jl")


# using .MiscFuncs, .AnalysisModule, .BBBModule, .Chi2Module, .KSModule
using .AnalysisModule

fName1 = "/home/shoram/Work/PhD_Thesis/Job17/Data/2vbb_Se82_Falaise_EneThetaPhi_1e8E.root"                           # root file for reference spectra
fName2 = "/home/shoram/Work/PhD_Thesis/Job17/Data/G0-G4_xi31_037_kappa_m06639_1e8E.root"   # root file for compared spectra

file1 = ROOTFile(fName1)
file2 = ROOTFile(fName2)

# singleElectronEnergies1 =
#     fill_from_root_file(file1, "tree", "reconstructedEnergy1") .+
#     fill_from_root_file(file1, "tree", "reconstructedEnergy2") # vector of single-electron energies for reference spectrum    
# singleElectronEnergies2 =
#     fill_from_root_file(file2, "tree", "reconstructedEnergy1") .+
#     fill_from_root_file(file2, "tree", "reconstructedEnergy2") # vector of single-electron energies for compared spectrum    

phi1 = fill_from_root_file(file1, "tree", "phi") # vector phi angles for reference spectrum    
phi2 = fill_from_root_file(file2, "tree", "phi") # vector phi angles for compared spectrum    

nSigma = 1
CL = 0.68
Δϕ = 15 # setting bin width
ΔE = 150

# sampleSizes = vcat(
#     collect(10_000:2000:18_000),
#     collect(20_000:10_000:100_000),
#     collect(100_000:50_000:400_000),
# )

# BBBPhi = BBB(phi1, phi2, 0, 180, Δϕ, nSigma)
# Chi2Phi = Chi2(phi1, phi2, 0, 180, Δϕ, sampleSizes, CL)
KSPhi = @suppress_err KS1(phi1, phi2, CL)


function KS1(
    vector1::Vector{T},
    vector2::Vector{T},
    CL::Real;
    nInARow = 3,
    nSamples = 100,
    sampleSizes = vcat(                 #default sample sizes, if different values are required, they can be passed as a new vector
        collect(10_000:1000:19_000),
        # collect(20_000:10_000:100_000),
        # collect(100_000:50_000:950_000),
    ),
    maxSampleSize = 100_000
) where T<:Real
    pVals = Vector{T}(undef, nSamples)   # initiaite a container to hold vectors of 100 p-values 
    efficiency = float(0)                # efficiency 
    nInARowCounter = 0                   # a counter for n in a row 100% efficiencies

    samples1 = Vector{T}(undef, nSamples)
    samples2 = Vector{T}(undef, nSamples)

    i = 1 
    sampleSize = sampleSizes[i]

    while (nInARowCounter != nInARow || sampleSize > maxSampleSize)
        if ( i < length(sampleSizes) )
            sampleSize = sampleSizes[i]
        else 
            sampleSize += 50_000 
        end
        samples1 = get_samples(vector1, sampleSize, nSamples, true)
        samples2 = get_samples(vector2, sampleSize, nSamples, true)

        @views pVals = ApproximateTwoSampleKSTest.(samples1, samples2) .|> pvalue
        efficiency = get_efficiency(pVals, CL, nSamples)

        println("efficiency: $efficiency")
        println("sampleSize: $sampleSize")

        if(efficiency == 1.0)
            nInARowCounter += 1
        else
            nInARowCounter = 0
        end
        i += 1
    end
    minEvents = sampleSize

    return KS(
        vector1,
        vector2,
        CL,
        minEvents
    )
end





BBBEne = BBB(singleElectronEnergies1, singleElectronEnergies2, 0, 3000, ΔE, nSigma)
Chi2Ene =
    Chi2(singleElectronEnergies1, singleElectronEnergies2, 450, 2100, ΔE, sampleSizes, CL)
KSEne = @suppress KS(
    singleElectronEnergies1,
    singleElectronEnergies2,
    0,
    3000,
    ΔE,
    sampleSizes,
    CL,
)

@show BBBPhi.minEvents
@show Chi2Phi.minEvents
@show KSPhi.minEvents

@show BBBEne.minEvents
@show Chi2Ene.minEvents
@show KSEne.minEvents
