


using MPO_MeanField
using CUDA
using ITensors, ITensorMPS
using Quantics, QuanticsTCI
import TensorCrossInterpolation as TCI
using LinearAlgebra
using Printf


println("CUDA Functional:", CUDA.functional())
if CUDA.functional()
    println("CUDA Device Name: ", CUDA.device())
    to_gpu = cu
    to_cpu = ITensors.cpu
    cleanup = CUDA.reclaim
else
    println("CUDA is not functional. Running on CPU.")
    to_gpu = identity
    to_cpu = identity
    cleanup = () -> nothing
end

println("="^50)
println("HartreeFockMPO — Timings in 2D")
println("="^50)


τ = (sqrt(10)-2.0)/2.0





function tx_qp(i)
     x,y = square_lattice_decoder(Int64(i),L)
     return -1.0 - V2 * (cos(2*pi * τ * (x-0.5)))
end
function ty_qp(i)
     x,y = square_lattice_decoder(Int64(i),L)
     return -1.0 - V2 * (cos(2*pi * τ * (y-0.5)))
end




function S(i)
     x,y = square_lattice_decoder(Int64(i),L)
     return W_Amp * (cos(pi*y) * cos(pi*x))
end


using ITensors
using CUDA
using HDF5
using Printf
using LinearAlgebra # For norm and tr

# Include your local modules here
# include("src/utils/quantics.jl")

function run_purification(L, V2, W_Amp, maxdim, tol)
    τ = (sqrt(10)-2.0)/2.0
    
    # 0-based decoding matching previous verification
    tx_qp(i) = -1.0 - V2 * (cos(2*pi * τ * (square_lattice_decoder(i-1, L)[1] - 0.5)))
    ty_qp(i) = -1.0 - V2 * (cos(2*pi * τ * (square_lattice_decoder(i-1, L)[2] - 0.5)))
    S(i) = W_Amp * (cos(pi * square_lattice_decoder(i-1, L)[2]) * cos(pi * square_lattice_decoder(i-1, L)[1]))

    params = ParametersSquare(L, (tx_qp, ty_qp), 0.3, nothing, S, 1e-6, tol, maxdim, 0.5, 200, 0.9, 0.1, 100)
    sys = System(params)

    t_rho0 = @elapsed ρ0 = construct_rho_0(sys, params, -3.0, 3.0)
    ρ0_gpu = to_gpu(ρ0)
    
    t_pur = @elapsed ρ_purified = perform_purification(ρ0_gpu, params; verbose=0)
    ρ_cpu = to_cpu(ρ_purified)
    
    return ρ_cpu, t_rho0, t_pur
end

function generate_table_1(L, V2, maxdim, tol)
    println("--- Table 1: MaxBond vs W_Amp (L=$L, V2=$V2) ---")
    
    W_amps = [1e-5, 0.1, 0.3, 0.5, 1.0]
    bond_dims = Float64[] # Float64 for homogeneous HDF5 matrix types
    pur_times = Float64[]
    
    for W in W_amps
        ρ_pur, _, t_pur = run_purification(L, V2, W, maxdim, tol)
        push!(bond_dims, maxlinkdim(ρ_pur))
        push!(pur_times, t_pur)
        
        @printf "W: %4.1f | MaxDim: %d | Time: %.2f s\n" W bond_dims[end] t_pur
        ρ_pur = nothing; GC.gc(); CUDA.functional() && CUDA.reclaim()
    end
    
    return hcat(W_amps, bond_dims, pur_times)
end

function generate_table_2(V2, W_Amp, maxdim, tol)
    println("\n--- Table 2: Time vs L (W_Amp=$W_Amp) ---")
    
    L_vals = [10, 12, 14, 16, 18, 20]
    sites  = Float64[]
    rho_times = Float64[]
    pur_times = Float64[]
    bond_dims = Float64[]
    
    for L in L_vals
        ρ_pur, t_rho0, t_pur = run_purification(L, V2, W_Amp, maxdim, tol)
        push!(sites, 2^L)
        push!(rho_times, t_rho0)
        push!(pur_times, t_pur)
        push!(bond_dims, maxlinkdim(ρ_pur))
        
        @printf "L: %d | Sites: %d | Rho0: %.2f s | Pur: %.2f s | MaxDim: %d\n" L (2^L) t_rho0 t_pur bond_dims[end]
        ρ_pur = nothing; GC.gc(); CUDA.functional() && CUDA.reclaim()
    end
    
    return hcat(L_vals, sites, rho_times, pur_times, bond_dims)
end

function generate_table_3(L, V2, W_Amp)
    println("\n--- Table 3: Error vs MaxDim (L=$L) ---")
    
    maxdims = [10, 20, 50, 100, 200,500,1000]
    tol = 1e-12 
    
    error_norms = Float64[]
    tr_vals = Float64[]
    
    # Placeholder for exact diagonalisation matrix
    # ρ_exact_matrix = solve_exact_diagonalization(...) 
    
    for md in maxdims
        ρ_pur, _, _ = run_purification(L, V2, W_Amp, md, tol)
        
        # NOTE: prod() scales exponentially with L. 
        # Only run this for small L (e.g., L <= 12).
        ρ_dense_itensor = prod(ρ_pur)
        ρ_mpo_matrix = matrix(ρ_dense_itensor)
        
        # error = norm(ρ_mpo_matrix - ρ_exact_matrix)
        error = 0.0 # Placeholder
        
        push!(error_norms, error)
        push!(tr_vals, real(tr(ρ_pur)))
        
        @printf "MaxDim: %d | Error: %.6e | Trace: %.6f\n" md error tr_vals[end]
        
        ρ_pur = nothing; ρ_dense_itensor = nothing; ρ_mpo_matrix = nothing
        GC.gc(); CUDA.functional() && CUDA.reclaim()
    end
    
    return hcat(maxdims, error_norms, tr_vals)
end

function main()
    # Parameters for the runs
    L_t1_t3 = 10
    V2      = 0.5
    W_Amp   = 0.3
    maxdim  = 500
    tol     = 1e-8
    
    # Generate data
    data_t1 = generate_table_1(L_t1_t3, V2, maxdim, tol)
    data_t2 = generate_table_2(V2, W_Amp, maxdim, tol)
    data_t3 = generate_table_3(L_t1_t3, V2, W_Amp)
    
    # Save to HDF5
    filename = "purification_results.h5"
    h5open(filename, "w") do file
        write(file, "table1_W_vs_MaxBond", data_t1)
        # Headers implicitly known: [W_amps, bond_dims, pur_times]
        
        write(file, "table2_L_vs_Time", data_t2)
        # Headers implicitly known: [L_vals, sites, rho_times, pur_times, bond_dims]
        
        write(file, "table3_MaxDim_vs_Error", data_t3)
        # Headers implicitly known: [maxdims, error_norms, tr_vals]
    end
    
    println("\nAll data successfully saved to $filename")
end

main()