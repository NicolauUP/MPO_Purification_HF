using MPO_MeanField

# Replace this single example with the parameter points for one scientific
# campaign. Do not generate a grid implicitly: every submitted point should be
# visible in version control.
campaign = (
    name="replace_with_campaign_name",
    runs=[
        (
            label="replace_with_physical_label",
            params=Parameters1D(
                L=2,
                t=-0.7,
                U=0.0,
                W=x -> (0.2, -0.1, 0.05, 0.4)[Int(x) + 1],
                S=nothing,
                tci_tol=1e-10,
                itensors_tol=1e-12,
                itensors_maxdim=64,
                density=0.5,
                purification_steps=35,
                scf_mixing=0.5,
                scf_tol=0.1,
                scf_max_iterations=10,
            ),
            spectral_bounds=(-5.0, 5.0),
            purification_method=:sp2,
            verify_spectral_bounds=true,
            verbose=:nothing,
        ),
    ],
)
