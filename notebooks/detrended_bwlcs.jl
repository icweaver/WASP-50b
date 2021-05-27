### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ b3ce86b4-bd60-11eb-17bf-ff4a35c62058
begin
	import PlutoUI as Pl
	using CairoMakie
	using Colors
	using CSV
	using DataFrames
	using Glob
	using Measurements
	using PyCall
	using Statistics, KernelDensity
	using Latexify, OrderedCollections
	using DataFramesMeta
	using NaturalSort
	using DelimitedFiles
end

# ╔═╡ 0158a760-1229-4089-bf90-7c7b2f1f548a
md"""
## Load data

First, let's load the relevant data needed for this notebook:
"""

# ╔═╡ 4b09c729-3395-4cee-bb69-bab59390845c
const DATA_DIR = "data/detrended/out_l_C/WASP50"

# ╔═╡ 100af59b-3a24-41d0-9cda-05592bd1778f
begin
	# Load IMACS
	cubes = Dict{String, Any}()
	
	for (i, (dirpath)) ∈ enumerate(sort(glob("$(DATA_DIR)/w50*IMACS*/wavelength")))
		fpaths = sort!(glob("$(dirpath)/wbin*/PCA_2/detrended_lc.dat"), lt=natural)
		
		dirpath_WLC = "$(dirname(dirpath))/white-light"
		
		# WLC BMA t₀
		t₀ = let
			df_BMA_WLC = CSV.File(
			"$(dirpath_WLC)/results.dat",
			comment = "#",
			normalizenames=true,
			) |> DataFrame
			df_BMA_WLC[findfirst(==("t0"), df_BMA_WLC.Variable), :].Value
		end
		
		# Pre-allocate matrices to hold detrended data and models in each bin
		N_time = size(readdlm("$(dirpath_WLC)/comps.dat", ','), 1)
		N_bins = length(fpaths)
		det_BLC_fluxes = Matrix{Float64}(undef, N_time, N_bins)
		det_BLC_models = copy(det_BLC_fluxes)
		
		# Populate each matrix
		for (fpath, lc, model) in zip(
			fpaths, eachcol(det_BLC_fluxes), eachcol(det_BLC_models)
			)
			df = CSV.File(
				fpath,
				header=["Time", "DetFlux", "DetFluxErr", "Model"],
				comment = "#",
				select=[:DetFlux, :Model],
			)
			lc .= df.DetFlux
			model .= df.Model
		end
				
		# Save
		cubes["Transit $i IMACS"] = Dict(
			"fluxes" => det_BLC_fluxes,
			"models" => det_BLC_models,
			"t₀" => t₀,
			"wbins" => readdlm("$(dirname(dirpath_WLC))/wbins.dat", comments=true),
			"tspec" => CSV.File(
				"$(dirname(dirpath_WLC))/transpec.csv",
				normalizenames = true,
			)
		)
	end
		
	cubes = sort(cubes)
end

# ╔═╡ ded314ba-1ebe-4f01-bd1b-652a0258f955
@bind transit Pl.Select(cubes.keys)

# ╔═╡ 2c0406e7-96e0-4a87-a91c-02d463e32ebc
md"""
## Helper functions
"""

# ╔═╡ f2da7123-cda9-47c3-aa72-4f47f4f8dfda
begin
	py"""
	import numpy as np
	import pickle
	
	def load_npz(fpath, allow_pickle=False):
		return np.load(fpath, allow_pickle=allow_pickle)[()]
	
	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	"""
	load_npz(s; allow_pickle=false) = py"load_npz"(s, allow_pickle=allow_pickle)
	load_pickle(s) = py"load_pickle"(s)
end;

# ╔═╡ 1dd4968e-959a-4f6e-a0e2-9fe9b8ecdd74
md"""
## Plot configs
"""

# ╔═╡ 0af97a94-cb08-40e2-8011-11c8696684fa
begin
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	const COLORS =  parse.(Colorant,
		[
			"#5daed9",  # Cyan
			"plum",
			"#f7ad4d",  # Yellow
			"mediumaquamarine",
			"#126399",  # Blue
			"#956cb4",  # Purple
			"#ff7f00",  # Orange
			"#029e73",  # Green
			"slategray",
		]
	)
end

# ╔═╡ bec88974-b150-4f53-9497-ddec4883ae17
function plot_BLCs(datas, models, wbins, errs; offset=0.3)
	median_prec = round(Int, median(errs))
	fig = Figure(resolution=FIG_TALL)
	ax_left = Axis(fig[1, 1], title = "Detrended BLCs")
	ax_right = Axis(fig[1, 2], title = "Residuals")
	ax_label = Axis(fig[1, 3], title = "Median precision: $(median_prec) ppm")
	axs = reshape(copy(fig.content), 1, 3)
	linkaxes!(axs...)
	
	# Color palette
	N_bins = size(datas, 2)
	resids = datas - models
	colors = to_colormap(:Spectral_4, N_bins) |> reverse
	
	# Arbitrary offsets for clarity
	offs = reshape(range(0, offset, length=N_bins), 1, :)
	baselines = ones(size(datas))
	for (data, model, resid, baseline, color, wbin, err) in zip(
			eachcol(datas .+ offs),
			eachcol(models .+ offs),
			eachcol(resids),
			eachcol(baselines .+ offs),
			colors,
			eachrow(wbins),
			errs,
		)	
		scatter!(ax_left, data, strokewidth=0, markersize=5, color=color)
		lines!(ax_left, model, linewidth=3, color=0.75*color)
		
		scatter!(ax_right, baseline + resid, markersize=5, color=color)
		lines!(ax_right, baseline, linewidth=3, color=0.75*color)
		text!(ax_label, "$(wbin[1]) - $(wbin[2]) Å, $(err) ppm";
			position = Point2f0(0, baseline[1]),
			textsize = 16,
			align = (:left, :center),
			offset = Point2f0(-10, 2),
			color = 0.75*color,
		)
	end
	
	hideydecorations!.(axs[:, 2:3], grid=false)
	hidespines!(axs[end])
	hidedecorations!(axs[end])
	#ylims!(ax_left, 0.95, 1.34)
	
	fig[1:2, 0] = Label(fig, "Relative flux + offset", rotation=π/2)
	fig[end, 2:3] = Label(fig, "Index")
	
	return fig
end

# ╔═╡ df1a160c-22ff-4c5e-a71f-b903d8a23ef1
let	
	tspec = cubes[transit]["tspec"]
	errs = maximum([tspec.Depthup_ppm_ tspec.DepthDown_ppm_], dims=2)
	
	p = plot_BLCs(
		cubes[transit]["fluxes"],
		cubes[transit]["models"],
		cubes[transit]["wbins"],
		round.(Int, errs),
	)
		
	p |> Pl.as_svg
end

# ╔═╡ Cell order:
# ╟─0158a760-1229-4089-bf90-7c7b2f1f548a
# ╠═4b09c729-3395-4cee-bb69-bab59390845c
# ╠═100af59b-3a24-41d0-9cda-05592bd1778f
# ╟─ded314ba-1ebe-4f01-bd1b-652a0258f955
# ╠═df1a160c-22ff-4c5e-a71f-b903d8a23ef1
# ╠═bec88974-b150-4f53-9497-ddec4883ae17
# ╟─2c0406e7-96e0-4a87-a91c-02d463e32ebc
# ╠═f2da7123-cda9-47c3-aa72-4f47f4f8dfda
# ╟─1dd4968e-959a-4f6e-a0e2-9fe9b8ecdd74
# ╠═0af97a94-cb08-40e2-8011-11c8696684fa
# ╠═b3ce86b4-bd60-11eb-17bf-ff4a35c62058
