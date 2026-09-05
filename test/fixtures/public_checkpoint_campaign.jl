using MPO_MeanField

campaign = CampaignSpec(
    name="checkpoint_smoke",
    cases=[CaseSpec(
        label="chain8",
        model=ChainModel(size=8, hopping=-0.7, interaction=0.0, filling=0.5),
        representation=QTTSettings(tci_tol=1e-10, cutoff=1e-12, maxdim=64),
        solver=SCFSettings(
            purification=:sp2,
            mixing=0.5,
            tolerance=0.1,
            maxiter=5,
            purification_maxiter=35,
            record_energy=false,
        ),
        spectral_bounds=(-3.0, 3.0),
    )],
)
