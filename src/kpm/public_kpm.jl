"""Sparse-vector KPM compatibility kernel for power-of-two lattice and graph models.

The polynomial recursion is independent of the lattice. The geometry-specific
part below only constructs the nearest-neighbour graph and maps local
estimators into Hartree/Fock fields. Physical storage is row-major with `x`
fastest; Hadamard colors use an interleaving of the coordinate bits so the
nested hierarchy remains two-dimensionally local for rectangles as well.
"""

struct KPMIterationRecord
    iteration::Int; spectral_lower::Float64; spectral_upper::Float64
    chemical_potential::Float64; trace::Float64; trace_error::Float64
    trace_squared_estimate::Float64; trace_idempotency_defect::Float64
    vh_residual::Float64; vf_residual::Float64; density_residual::Float64
    bond_residual::Float64; two_cycle_residual::Float64; energy_total::Float64
    checkerboard_order::Float64
end
struct KPMDiagnostics
    history::Vector{KPMIterationRecord}; converged::Bool; termination_reason::Symbol; audit::NamedTuple
end
_kpm_change(a, b) = norm(a - b) / max(norm(b), sqrt(eps(Float64)))

function _kpm_coefficients(M::Int, μ::Float64)
    -1 < μ < 1 || throw(ArgumentError("scaled KPM chemical potential must lie in (-1,1)"))
    θ, angle = acos(μ), pi / (M + 1)
    c = Vector{Float64}(undef, M + 1); c[1] = 2 * (pi - θ) / pi
    for m in 1:M
        c[m + 1] = -2sin(m * θ) / (pi * m)
    end
    return [c[m + 1] * ((M - m + 1) * cos(m * angle) + sin(m * angle) * cot(angle)) / (M + 1) for m in 0:M]
end

function _kpm_coordinate_code(x::Int, y::Int, x_bits::Int, y_bits::Int)
    code = 0
    position = 0
    for bit in 0:(min(x_bits, y_bits) - 1)
        code |= ((y >> bit) & 1) << position
        position += 1
        code |= ((x >> bit) & 1) << position
        position += 1
    end
    for bit in min(x_bits, y_bits):(x_bits - 1)
        code |= ((x >> bit) & 1) << position
        position += 1
    end
    for bit in min(x_bits, y_bits):(y_bits - 1)
        code |= ((y >> bit) & 1) << position
        position += 1
    end
    return code
end

function _kpm_hadamard(codes::AbstractVector{<:Integer}, R::Int, seed::Int)
    N = length(codes)
    ispow2(N) && 0 < R <= N || throw(ArgumentError("Hadamard KPM probes require 0 < R ≤ N and power-of-two N"))
    sort!(collect(Int, codes)) == collect(0:(N - 1)) || throw(ArgumentError(
        "Hadamard row codes must be a zero-based permutation of 0:N-1",
    ))
    signs = ifelse.(rand(Xoshiro(seed), Bool, N), 1.0, -1.0)
    probes = Matrix{Float64}(undef, N, R)
    for column in 0:(R - 1), row in eachindex(codes)
        probes[row, column + 1] = isodd(count_ones(codes[row] & xor(column, column >> 1))) ? -signs[row] : signs[row]
    end
    probes
end

function _kpm_apply(H, Z, c; synchronize=() -> nothing)
    previous, current, following, result = copy(Z), similar(Z), similar(Z), similar(Z)
    mul!(current, H, Z); @. result = (c[1] / 2) * previous + c[2] * current
    for m in 2:(length(c) - 1)
        mul!(following, H, current); @. following = 2 * following - previous; @. result += c[m + 1] * following
        previous, current, following = current, following, previous
    end
    synchronize(); result
end

function _kpm_moments(H, Z, M; synchronize=() -> nothing)
    moment = Vector{Float64}(undef, M + 1); normalization = inv(size(Z, 2))
    previous, current, following = copy(Z), similar(Z), similar(Z)
    moment[1] = real(dot(vec(Z), vec(previous))) * normalization
    mul!(current, H, Z); moment[2] = real(dot(vec(Z), vec(current))) * normalization
    for m in 2:M
        mul!(following, H, current); @. following = 2 * following - previous
        moment[m + 1] = real(dot(vec(Z), vec(following))) * normalization
        previous, current, following = current, following, previous
    end
    synchronize(); moment
end

function _kpm_mu(moments, Ne; tolerance)
    trace(μ) = begin c = _kpm_coefficients(length(moments) - 1, μ); c[1] * moments[1] / 2 + dot(@view(c[2:end]), @view(moments[2:end])) end
    lower, upper = -1.0 + 1e-10, 1.0 - 1e-10
    trace(lower) <= Ne <= trace(upper) || throw(ArgumentError("KPM trace bracket does not contain target filling"))
    μ, t = 0.0, 0.0
    for _ in 1:80
        μ = (lower + upper) / 2; t = trace(μ); abs(t - Ne) <= tolerance && break
        t < Ne ? (lower = μ) : (upper = μ)
    end
    (; scaled_mu=μ, trace=t)
end

function _kpm_data(model::SquareModel)
    nx, ny = model.size
    N = nx * ny
    x_bits, y_bits = trailing_zeros(nx), trailing_zeros(ny)
    onsite, seed, bonds, hopping = zeros(N), zeros(N), Tuple{Int,Int,Symbol}[], Float64[]
    codes = Vector{Int}(undef, N)
    site_index(x, y) = x + nx * y + 1
    tx(x,y) = model.hopping[1] isa Number ? Float64(model.hopping[1]) : Float64(model.hopping[1](x,y))
    ty(x,y) = model.hopping[2] isa Number ? Float64(model.hopping[2]) : Float64(model.hopping[2](x,y))
    for y in 0:(ny - 1), x in 0:(nx - 1)
        i = site_index(x, y)
        onsite[i] = isnothing(model.potential) ? 0.0 : Float64(model.potential(x, y))
        seed[i] = isnothing(model.seed) ? 0.0 : Float64(model.seed(x, y))
        codes[i] = _kpm_coordinate_code(x, y, x_bits, y_bits)
        if x < nx - 1; push!(bonds, (i, site_index(x + 1, y), :horizontal)); push!(hopping, tx(x,y)); end
        if y < ny - 1; push!(bonds, (i, site_index(x, y + 1), :vertical)); push!(hopping, ty(x,y)); end
    end
    (; nx, ny, x_bits, y_bits, N, onsite, seed, bonds, hopping, codes)
end

function _kpm_data(model::GraphModel)
    N = size(model.hopping, 1)
    bonds = Tuple{Int,Int,Symbol}[]
    hopping = Float64[]
    for column in 1:N
        for pointer in nzrange(model.hopping, column)
            row = model.hopping.rowval[pointer]
            row < column || continue
            push!(bonds, (row, column, :graph))
            push!(hopping, model.hopping.nzval[pointer])
        end
    end
    return (; nx=nothing, ny=nothing, x_bits=nothing, y_bits=nothing, N,
        onsite=copy(model.potential), seed=copy(model.seed), bonds, hopping,
        codes=copy(model.probe_codes))
end

function _kpm_hamiltonian(data, VH, VF)
    rows, columns, values = Int[], Int[], Float64[]
    for i in 1:data.N; push!(rows,i); push!(columns,i); push!(values,data.onsite[i]+VH[i]); end
    for k in eachindex(data.bonds)
        i,j,_ = data.bonds[k]; v=data.hopping[k]+VF[k]; append!(rows,(i,j)); append!(columns,(j,i)); append!(values,(v,v))
    end
    sparse(rows, columns, values, data.N, data.N)
end
function _kpm_bounds(data, VH, VF)
    radius=zeros(data.N)
    for k in eachindex(data.bonds); i,j,_=data.bonds[k]; v=abs(data.hopping[k]+VF[k]); radius[i]+=v; radius[j]+=v; end
    diagonal=data.onsite+VH; minimum(diagonal-radius)-0.1, maximum(diagonal+radius)+0.1
end
function _kpm_local(PZ,Z,bonds)
    R=size(Z,2); density=vec(sum(PZ .* Z;dims=2))./R; order=Vector{Float64}(undef,length(bonds))
    for k in eachindex(bonds); i,j,_=bonds[k]; order[k]=(dot(@view(PZ[i,:]),@view(Z[j,:]))+dot(@view(PZ[j,:]),@view(Z[i,:])))/(2R); end
    tr=sum(PZ .* Z)/R; tr2=sum(abs2,PZ)/R
    (; density, order, trace_squared=tr2, idempotency=abs(tr-tr2)/max(abs(tr),sqrt(eps(Float64))))
end
function _kpm_fields(density, order, data, U)
    VH=zeros(data.N); for (i,j,_) in data.bonds; VH[i]+=U*density[j]; VH[j]+=U*density[i]; end; VH, -U .* order
end
function _kpm_energy(data,density,order,U)
    kinetic=dot(data.onsite,density)+2dot(data.hopping,order); hartree=sum(U*density[i]*density[j] for (i,j,_) in data.bonds); fock=-sum(U*abs2(v) for v in order)
    (; kinetic, hartree, fock, interaction=hartree+fock, total=kinetic+hartree+fock)
end
function _kpm_observables(data,density,order,energy,idem)
    hb=Tuple{Int,Int}[]; vb=Tuple{Int,Int}[]; ho=ComplexF64[]; vo=ComplexF64[]
    for (bond,v) in zip(data.bonds,order)
        if bond[3] == :horizontal; push!(hb,(bond[1],bond[2]));push!(ho,v)
        elseif bond[3] == :vertical; push!(vb,(bond[1],bond[2]));push!(vo,v); end
    end
    (; site_density=density, bonds=data.bonds, bond_order=ComplexF64.(order),
        horizontal_bonds=hb,vertical_bonds=vb,horizontal_bond_order=ho,vertical_bond_order=vo,
        particle_number=sum(density),energy,hermiticity_residual=NaN,
        idempotency_residual=idem,stationarity_residual=NaN)
end

function _solve_kpm(model::Union{SquareModel,GraphModel}, solver::SCFSettings, settings::KPMSettings, runtime::RuntimeReport)
    data=_kpm_data(model); Ne=round(Int,model.filling*data.N); cuda=nothing; synchronize=()->nothing
    if runtime.active_backend == :cuda; cuda,error=_cuda_module_or_error(); isnothing(error)||throw(ArgumentError(error)); synchronize=()->_cuda_call(cuda,:synchronize); end
    hostZ=_kpm_hadamard(data.codes,settings.probes,settings.probe_seed); Z=isnothing(cuda) ? hostZ : _cuda_call(cuda,:CuArray,hostZ)
    VH=copy(data.seed); VF=zeros(length(data.bonds)); previous_density=nothing; previous_order=nothing; two_steps=nothing; history=KPMIterationRecord[]; stable=0; termination=:max_iterations; last=nothing
    trace_tolerance=max(1e-6*Ne,1e-6); tolerance=solver.tolerance/100
    for iter in 1:solver.maxiter
        lo,hi=_kpm_bounds(data,VH,VF); center=(lo+hi)/2; half=(hi-lo)/2; scaled=(_kpm_hamiltonian(data,VH,VF)-center*sparse(I,data.N,data.N))/half
        H=isnothing(cuda) ? scaled : Base.invokelatest(getproperty(getproperty(cuda,:CUSPARSE),:CuSparseMatrixCSR),scaled)
        μ = _kpm_mu(_kpm_moments(H, Z, settings.moments; synchronize=synchronize), Ne; tolerance=trace_tolerance)
        PZ = _kpm_apply(H, Z, _kpm_coefficients(settings.moments, μ.scaled_mu); synchronize=synchronize)
        estimates = _kpm_local(isnothing(cuda) ? PZ : Array(PZ), hostZ, data.bonds)
        newVH,newVF=_kpm_fields(estimates.density,estimates.order,data,Float64(model.interaction)); vh=isnothing(previous_density) ? Inf : _kpm_change(newVH,VH); vf=isnothing(previous_density) ? Inf : _kpm_change(newVF,VF); dr=isnothing(previous_density) ? Inf : _kpm_change(estimates.density,previous_density); br=isnothing(previous_order) ? Inf : _kpm_change(estimates.order,previous_order); cycle=isnothing(two_steps) ? Inf : _kpm_change(estimates.density,two_steps); energy=_kpm_energy(data,estimates.density,estimates.order,Float64(model.interaction)); checker=model isa SquareModel ? sum((iseven(x+y) ? 1.0 : -1.0)*estimates.density[x + data.nx*y + 1] for y in 0:(data.ny-1),x in 0:(data.nx-1))/data.N : NaN
        push!(history,KPMIterationRecord(iter,lo,hi,center+half*μ.scaled_mu,sum(estimates.density),abs(sum(estimates.density)-Ne),estimates.trace_squared,estimates.idempotency,vh,vf,dr,br,cycle,energy.total,checker)); last=(;density=estimates.density,order=estimates.order,energy,idem=estimates.idempotency,stationarity=max(vh,vf,dr,br)); stable=all(isfinite,(vh,vf,dr,br))&&max(vh,vf,dr,br)<=tolerance ? stable+1 : 0
        if stable >= 2
            termination = :converged
            break
        elseif iter >= 3 && cycle <= tolerance && dr > tolerance
            termination = :two_cycle_detected
            break
        end
        two_steps=previous_density;previous_density=estimates.density;previous_order=estimates.order; if iter==1;VH,VF=newVH,newVF else;VH=solver.mixing.*newVH.+(1-solver.mixing).*VH;VF=solver.mixing.*newVF.+(1-solver.mixing).*VF;end
    end
    auditZhost = _kpm_hadamard(data.codes, settings.audit_probes, settings.audit_seed)
    auditZ = isnothing(cuda) ? auditZhost : _cuda_call(cuda, :CuArray, auditZhost)
    lo, hi = _kpm_bounds(data, VH, VF); center = (lo + hi) / 2; half = (hi - lo) / 2
    scaled = (_kpm_hamiltonian(data, VH, VF) - center * sparse(I, data.N, data.N)) / half
    H = isnothing(cuda) ? scaled : Base.invokelatest(getproperty(getproperty(cuda, :CUSPARSE), :CuSparseMatrixCSR), scaled)
    μ = _kpm_mu(_kpm_moments(H, auditZ, settings.audit_moments; synchronize=synchronize), Ne; tolerance=trace_tolerance)
    auditPZ = _kpm_apply(H, auditZ, _kpm_coefficients(settings.audit_moments, μ.scaled_mu); synchronize=synchronize)
    auditlocal = _kpm_local(isnothing(cuda) ? auditPZ : Array(auditPZ), auditZhost, data.bonds)
    audit=(;moments=settings.audit_moments,probes=settings.audit_probes,seed=settings.audit_seed,trace=sum(auditlocal.density),trace_error=abs(sum(auditlocal.density)-Ne),trace_idempotency_defect=auditlocal.idempotency,density_max_abs_difference=maximum(abs,auditlocal.density-last.density),density_rms_difference=norm(auditlocal.density-last.density)/sqrt(data.N),bond_max_abs_difference=maximum(abs,auditlocal.order-last.order),bond_rms_difference=norm(auditlocal.order-last.order)/sqrt(length(data.bonds)))
    diagnostics=KPMDiagnostics(history,termination==:converged,termination,audit); (;converged=diagnostics.converged,termination_reason=termination,diagnostics,observables=_kpm_observables(data,last.density,last.order,last.energy,last.idem))
end
