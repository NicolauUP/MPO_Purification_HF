using MPO_MeanField

# Canonical public campaign template. Replace the example values with one
# explicit, scientifically justified case at a time. Function-valued hopping,
# potential, and seed terms may be defined above the CampaignSpec.
campaign = CampaignSpec(
    name="replace_with_campaign_name",
    cases=[
        CaseSpec(
            label="replace_with_physical_label",
            model=ChainModel(
                size=4,
                hopping=-0.7,
                interaction=0.0,
                potential=x -> (0.2, -0.1, 0.05, 0.4)[Int(x) + 1],
                seed=nothing,
                filling=0.5,
            ),
            representation=QTTSettings(
                encoding=:binary,
                tci_tol=1e-10,
                cutoff=1e-12,
                maxdim=64,
            ),
            solver=SCFSettings(
                purification=:sp2,
                mixing=0.5,
                tolerance=0.1,
                maxiter=10,
                purification_maxiter=35,
            ),
            spectral_bounds=(-5.0, 5.0),
        ),
    ],
)
