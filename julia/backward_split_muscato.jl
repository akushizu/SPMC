using Pkg
Pkg.activate(".")
Pkg.add("JLD2")
Pkg.add("CairoMakie")
Pkg.add("FFTW")

using JLD2
using SparseArrays
using LinearAlgebra
using FFTW
using CairoMakie

function backward()

    hb = 6.62607015 / 1.602176634 / 2 / pi
    m = 0.32
    e = 1.0
    k0 = 0.7
    s0 = 2.852
    x0 = -15.0
    xc = 0.0
    a = 0.3
    s = 1.0
    
    # -------- Step size ----------------------
    dx = 0.005
    dt = 0.001
    l = (0.5 * im * hb / m) * (dt / dx^2)
    
    x = collect(-40 : dx : 40)
    t = collect(0 : dt : 25)
    
    # ------ Initial & boundary condition ------------------
    w = sqrt(sqrt(1 / (2 * pi * s0^2))) * exp.(-(x .- x0).^2 / (4 * s0^2)) .* exp.(im * k0 * (x .- xc))
    
    N = 1 #sqrt( dx * sum( abs(w).^2 ) )

    w = w / N
    
    # --------- Potential ------------------------
    Vw = e * a * exp.(-(x .- xc).^2 / (2 * s^2))
    
    # --------- Backward difference ------------------------
    W = spdiagm(0 => (1 + 2 * l) .+ (im * dt / hb) .* Vw, 1 => fill(-l, length(x) - 1), -1 => fill(-l, length(x) - 1))
    
    DataType = NamedTuple{(:t, :x, :Amp, :time), Tuple{Float64, Vector{Float64}, Vector{ComplexF64}, Float64}}
    adata = Vector{DataType}(undef, Int(((length(t) - 1) / 50) + 1))
    adata[1] = (t = t[1], x = x, Amp = deepcopy(w), time = 0.)
    
    t0 = time_ns()

    for i in eachindex(t)[2:end]

        w = W \ w

        N = 1 #sqrt( dx * sum( abs(w1).^2 ) )
        
        w .= w ./ N
        
        if mod(i - 1, 50) == 0

            adata[Int(( i - 1 ) / 50 + 1)] = (t = t[i], x = x, Amp = deepcopy(w), time = time() - t0)
            
        end

    end

	println("Elapsed time is $((time_ns()-t0)/1e9) seconds.")

    # ----------save----------------------------
    jldsave("adata.jld2"; adata)
    
    # ---------- Plotting & Video ----------------------------
    colors = Makie.wong_colors()
    fig = Figure()
    
    amp = Observable(adata[1].Amp)
    tt = Observable(adata[1].t)
    
    ax1 = Axis(fig[1, 1],
            title = @lift("den, t = $(round($tt, digits=2)) fs"),
            xlabel="x")
    lines!(ax1, x, 0.5 .* Vw, linestyle=:dash, color=colors[2])
    lines!(ax1, x, @lift(abs2.($amp)), color=colors[1])
    xlims!(ax1, x[1], x[end])

    ax2 = Axis(fig[1, 2],
            title="abs",
            xlabel="x")
    lines!(ax2, x, Vw, linestyle=:dash, color=colors[2])
    lines!(ax2, x, @lift(abs.($amp)), color=colors[1])
    xlims!(ax2, x[1], x[end])

    ax3 = Axis(fig[2, 1],
            title="real",
            xlabel="x")
    lines!(ax3, x, Vw, linestyle=:dash, color=colors[2])
    lines!(ax3, x, @lift(abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax3, x, @lift(-abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax3, x, @lift(real.($amp)), color=colors[1])
    xlims!(ax3, x[1], x[end])

    ax4 = Axis(fig[2, 2],
            title="imag",
            xlabel="x")
    lines!(ax4, x, Vw, linestyle=:dash, color=colors[2])
    lines!(ax4, x, @lift(abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax4, x, @lift(-abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax4, x, @lift(imag.($amp)), color=colors[1])
    xlims!(ax4, x[1], x[end])
    
    record(fig, "wob.mp4", 1:length(adata); framerate=12) do i
        
        amp[] = adata[i].Amp
        tt[] = adata[i].t
        autolimits!(ax1)
        autolimits!(ax2)
        autolimits!(ax3)
        autolimits!(ax4)

    end

    return
end

function split_op()

    hb = 6.62607015 / 1.602176634 / 2 / pi
    m = 0.32
    e = 1.0
    k0 = 0.7
    s0 = 2.852
    x0 = -15.0
    xc = 0.0
    a = 0.3
    s = 1.0
    
    # -------- Step size ----------------------
    N = 16384
    dt = 0.001
    
    x = collect(LinRange(-40, 40, N))
    dx = x[2] - x[1]
    k = collect(LinRange(-pi/dx, pi/dx, N))
    t = collect(0 : dt : 25)
    
    # ------ Initial & boundary condition ------------------
    psi = sqrt(sqrt(1 / (2 * pi * s0^2))) * exp.(-(x .- x0).^2 / (4 * s0^2)) .* exp.(im * k0 * (x .- xc))
    V = e * a * exp.(-(x .- xc).^2 / (2 * s^2))
    
    Ut = exp.(-im * (hb^2 * k.^2 / m / 2) * dt / hb)
    Uv = exp.(-im * V * dt / 2 / hb)
    
    DataType = NamedTuple{(:t, :x, :Amp, :time), Tuple{Float64, Vector{Float64}, Vector{ComplexF64}, Float64}}
    edata = Vector{DataType}(undef, Int(((length(t) - 1) / 50) + 1))
    edata[1] = (t = t[1], x = x, Amp = deepcopy(psi), time = 0.)
    
    p = plan_fft(psi)
    ip = plan_ifft(psi)
    
    t0 = time_ns()

    for i in eachindex(t)[2:end]

        psi .= Uv .* psi
        psih = fftshift(p * psi)
        psih .= Ut .* psih
        psi .= ip * ifftshift(psih)
        psi .= Uv .* psi
        
        if mod(i - 1, 50) == 0

            edata[Int(( i - 1 ) / 50 + 1)] = (t = t[i], x = x, Amp = deepcopy(psi), time = time() - t0)

        end

    end

    println("Elapsed time is $((time_ns()-t0)/1e9) seconds.")

    # ----------save----------------------------
    jldsave("edata.jld2"; edata)
    
    # ---------- Plotting & Video ----------------------------
    colors = Makie.wong_colors()
    fig = Figure()
    
    amp = Observable(edata[1].Amp)
    tt = Observable(edata[1].t)
    
    ax1 = Axis(fig[1, 1],
            title = @lift("den, t = $(round($tt, digits=2)) fs"),
            xlabel="x")
    lines!(ax1, x, 0.5 .* V, linestyle=:dash, color=colors[2])
    lines!(ax1, x, @lift(abs2.($amp)), color=colors[1])
    xlims!(ax1, x[1], x[end])

    ax2 = Axis(fig[1, 2],
            title="abs",
            xlabel="x")
    lines!(ax2, x, V, linestyle=:dash, color=colors[2])
    lines!(ax2, x, @lift(abs.($amp)), color=colors[1])
    xlims!(ax2, x[1], x[end])

    ax3 = Axis(fig[2, 1],
            title="real",
            xlabel="x")
    lines!(ax3, x, V, linestyle=:dash, color=colors[2])
    lines!(ax3, x, @lift(abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax3, x, @lift(-abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax3, x, @lift(real.($amp)), color=colors[1])
    xlims!(ax3, x[1], x[end])

    ax4 = Axis(fig[2, 2],
            title="imag",
            xlabel="x")
    lines!(ax4, x, V, linestyle=:dash, color=colors[2])
    lines!(ax4, x, @lift(abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax4, x, @lift(-abs.($amp)), color=:black, linestyle=:dot)
    lines!(ax4, x, @lift(imag.($amp)), color=colors[1])
    xlims!(ax4, x[1], x[end])
    
    record(fig, "wos.mp4", 1:length(edata); framerate=12) do i

        amp[] = edata[i].Amp
        tt[] = edata[i].t
        autolimits!(ax1)
        autolimits!(ax2)
        autolimits!(ax3)
        autolimits!(ax4)
        
    end

    return
end

function compare()

    adata = load("adata.jld2", "adata")
    edata = load("edata.jld2", "edata")

    e = 1.0
    xc = 0.0
    a = 0.12
    s = 1.0
    
    x = collect(-40:0.005:40)
    V = e * a * exp.(-(x .- xc).^2 / (2 * s^2))
    
    colors = Makie.wong_colors()
    fig = Figure()
    
    aAmp = Observable(adata[1].Amp)
    eAmp = Observable(edata[1].Amp)
    tt = Observable(adata[1].t)
    
    ax = Axis(fig[1, 1], title = @lift("den, t = $(round($tt, digits=2)) fs"), xlabel="x")
    
    lines!(ax, x, 0.5 .* V, linestyle=:dash, color=colors[3])
    l1 = lines!(ax, adata[1].x, @lift(abs2.($aAmp)), color=colors[1], label="Backward difference")
    l2 = lines!(ax, edata[1].x, @lift(abs2.($eAmp)), color=colors[2], label="Split operator")
    
    xlims!(ax, x[1], x[end])
    axislegend(ax)
    
    record(fig, "woc.mp4", 1:length(adata); framerate=12) do i

        aAmp[] = adata[i].Amp
        eAmp[] = edata[i].Amp
        tt[] = adata[i].t
        autolimits!(ax)
        
    end
end

function main()
    
    backward()
    
    split_op()
    
    compare()
    
end

main()