### THIS SCRIPT IS USED TO CALCULATE PVALUES FOR AD TEST. WATCH OUT IT TAKES VERY LONG TIME TO FINISH!!

using Pkg;
Pkg.activate(Base.current_project());

using Revise
using StatsPlots, UnROOT, DataFrames, CSV

push!(LOAD_PATH, pwd() * "/src/")
using AnalysisModule
Revise.track(AnalysisModule)

## load data
const fName1 = "../Job17/Data/2vbb_Se82_Falaise_EneThetaPhi_1e8E.root"       # root file for reference spectra
const fName2 = "../Job17/Data/G0-G4_xi31_037_kappa_m06639_1e8E.root"         # root file for compared spectra

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
const CLs = [0.68, 0.90, 0.95]
const sampleSizes = vcat(collect(20_000:10_000:100_000), collect(150_000:50_000:1_000_000))

function main()
    for CL in CLs
        @show "First CL starting, CL = $CL"

        
        Chi2Phi = Chi2( phi1, phi2, CL, 0, 180, 15, sampleSizes = sampleSizes )
        Chi2Ene = Chi2( singleElectronEnergies1, singleElectronEnergies2, CL, 450, 2100, 150, sampleSizes = sampleSizes )
        @show Chi2Phi.minEvents, Chi2Ene.minEvents

        KSPhi = KS( phi1, phi2, CL, sampleSizes = sampleSizes )
        KSEne = KS( singleElectronEnergies1, singleElectronEnergies2, CL, sampleSizes = sampleSizes )
        @show KSPhi.minEvents, KSEne.minEvents

        ADPhi = AD( phi1, phi2, CL, sampleSizes = sampleSizes )
        ADEne = AD( singleElectronEnergies1, singleElectronEnergies2, CL, sampleSizes = sampleSizes )
        @show ADPhi.minEvents, ADEne.minEvents

        @show "Starting get_pVals, CL = $CL" 
        pValsChi2Phi = get_pVals(Chi2Phi, sampleSizes)
        pValsChi2Ene = get_pVals(Chi2Ene, sampleSizes)
        pValsKSPhi = get_pVals(KSPhi, sampleSizes)
        pValsKSEne = get_pVals(KSEne, sampleSizes)
        pValsADPhi = get_pVals(ADPhi, sampleSizes)
        pValsADEne = get_pVals(ADEne, sampleSizes)
        @show "Finished get_pVals, CL = $CL" 


        dfChi2Phi = DataFrame(pValsChi2Phi, ["$sampleSize" for sampleSize in sampleSizes])
        dfChi2Ene = DataFrame(pValsChi2Ene, ["$sampleSize" for sampleSize in sampleSizes])    
        dfKSPhi = DataFrame(pValsKSPhi, ["$sampleSize" for sampleSize in sampleSizes])
        dfKSEne = DataFrame(pValsKSEne, ["$sampleSize" for sampleSize in sampleSizes])
        dfADPhi = DataFrame(pValsADPhi, ["$sampleSize" for sampleSize in sampleSizes])
        dfADEne = DataFrame(pValsADEne, ["$sampleSize" for sampleSize in sampleSizes])

        saveNameChi2Phi = string("pValseChi2Phi_CL=$CL.csv")
        saveNameChi2Ene = string("pValseChi2Ene_CL=$CL.csv")
        saveNameKSPhi = string("pValseKSPhi_CL=$CL.csv")
        saveNameKSEne = string("pValseKSEne_CL=$CL.csv")
        saveNameADPhi = string("pValseADPhi_CL=$CL.csv")
        saveNameADEne = string("pValseADEne_CL=$CL.csv")

        CSV.write(saveNameChi2Phi, dfChi2Phi)
        CSV.write(saveNameChi2Ene, dfChi2Ene)    
        CSV.write(saveNameKSPhi, dfKSPhi)
        CSV.write(saveNameKSEne, dfKSEne)
        CSV.write(saveNameADPhi, dfADPhi)
        CSV.write(saveNameADEne, dfADEne)
    end
end


main()