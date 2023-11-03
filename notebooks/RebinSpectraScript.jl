using DataFramesMeta, CSV, DelimitedFiles

push!(LOAD_PATH, "/home/shoram/.julia/dev/MPRebinSpectra/") # add package to path so it can be used by using MPRebinSpectra
push!(LOAD_PATH, "/home/shoram/Work/PhD_Thesis/SNImproved2bbSpectra/src") # add package to path so it can be used by using MPRebinSpectra
using AnalysisModule
using MPRebinSpectra

begin
    const G0 = 0.334863e-46 # [MeV]
    const G2 = 0.148350E-46 # [MeV]
    const G22 = 0.188606E-47 # [MeV]
    const G4 = 0.824467E-47 # [MeV]

    const H0 = 0.226022E-46 # [MeV]
    const H2 = 0.929409E-47 # [MeV]
    const H22 = 0.108907E-47 # [MeV]
    const H4 = 0.484671E-47 # [MeV]

end

ξ51 = 0.1397 # for SSD
ξ31 = 0.3738 # declare xi31
Kappa = round(get_kappa(ξ31), digits = 4)

begin
    dataDir = "/home/shoram/Work/PhD_Thesis/Data_RebinnedSpectra/Se82/2vbb_Angular/FourOrders/";

    G0file = "/home/shoram/Work/PhD_Thesis/Data_RebinnedSpectra/Se82/1-Gfactors/spectrumG0.dat"
    G2file = "/home/shoram/Work/PhD_Thesis/Data_RebinnedSpectra/Se82/1-Gfactors/spectrumG2.dat"
    G22file = "/home/shoram/Work/PhD_Thesis/Data_RebinnedSpectra/Se82/1-Gfactors/spectrumG22.dat"
    G4file = "/home/shoram/Work/PhD_Thesis/Data_RebinnedSpectra/Se82/1-Gfactors/spectrumG4.dat"
end

#######################################################
########## FOR COMBINING THE INDIVIDUAL Gi FILES!!!!!!!
#######################################################

function add_spectra(ξ31, ξ51, G0file, G2file, G22file, G4file)
    dG0 = readdlm(G0file, Float64)[:, 3]
    dG2 = readdlm(G2file, Float64)[:, 3]
    dG22 = readdlm(G22file, Float64)[:, 3]
    dG4 = readdlm(G4file, Float64)[:, 3]

    dG = @. dG0 + ξ31 * dG2 + 1 / 3 * (ξ31)^2 * dG22 + (1 / 3 * (ξ31)^2 + ξ51) * dG4
    return dG
end


dG = add_spectra(ξ31, ξ51, G0file, G2file, G22file, G4file);

let E1 = readdlm(G0file, Float64)[:, 1], E2 = readdlm(G0file, Float64)[:, 2]
    saveFileName = "xi51=$(ξ51)_xi31=$(ξ31)_K2v=$(Kappa)_G0_G2_G22_G4.csv";
    
    open(joinpath(dataDir, saveFileName), "w") do io
        writedlm(io, [E1 E2 dG], "    ")
    end

    df_raw =
    CSV.File(
        joinpath(dataDir, saveFileName),
        delim = "    ",
        header = ["E1", "E2", "dGdE"],
    ) |> DataFrame

    df = rebin2D(df_raw, 0.001) |> normalize2D! |> get_cdf!
    CSV.write(joinpath(dataDir, "rebinned", saveFileName), df)
end;

#######################################################
########## FOR WORKING WITH FULL DATASET!!!!!!!!!!!!!!!
#######################################################
let E1 = readdlm(G0file, Float64)[:, 1], E2 = readdlm(G0file, Float64)[:, 2]
    inFileName = "FromRasto_xi31=0.3738_K=-0.6638.dat"
    saveFileName = "xi51=$(ξ51)_xi31=$(ξ31)_K2v=$(Kappa)_fromFullDataSet.csv";

    df_raw =
        CSV.File(
            joinpath(dataDir, inFileName),
            delim = "    ",
            header = ["E1", "E2", "dGdE"],
        ) |> DataFrame

    df = rebin2D(df_raw, 0.001) |> normalize2D! |> get_cdf!
    CSV.write(joinpath(dataDir, "rebinned", saveFileName), df)
end;