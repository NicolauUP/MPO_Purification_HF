@testset "M1.1 basis and matrix conventions" begin
    sites = siteinds("Qubit", 2)
    bra_states, ket_states = MPO_MeanField.precompute_qtt_states(sites)

    matrix_element(mpo, i, j) = MatrixChecker(
        mpo, sites, i, j, bra_states, ket_states,
    )
    matrix_of(mpo) = [matrix_element(mpo, i, j) for i in 1:4, j in 1:4]

    # MatrixChecker uses one-based Julia indices for the zero-based binary
    # states: entry (i,j) is <i-1|M|j-1>.
    first_bit = OpSum()
    first_bit += 1.0, "P+", 1
    second_bit = OpSum()
    second_bit += 1.0, "P+", 2
    @test matrix_of(MPO(first_bit, sites)) == Diagonal([0.0, 0.0, 1.0, 1.0])
    @test matrix_of(MPO(second_bit, sites)) == Diagonal([0.0, 1.0, 0.0, 1.0])

    T_R, T_L = MPO_MeanField.build_translation_chain(sites)
    expected_T_R = [
        0.0 1.0 0.0 0.0
        0.0 0.0 1.0 0.0
        0.0 0.0 0.0 1.0
        0.0 0.0 0.0 0.0
    ]
    expected_T_L = transpose(expected_T_R)

    # Full-matrix equality checks every element, including both open-chain
    # corners where a periodic wraparound term would occur.
    @test matrix_of(T_R) == expected_T_R
    @test matrix_of(T_L) == expected_T_L
    @test matrix_of(T_L) == matrix_of(T_R)'
end
