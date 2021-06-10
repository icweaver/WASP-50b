### A Pluto.jl notebook ###
# v0.14.8

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

# ╔═╡ 818783f8-7164-466e-b5a7-b75eaefe6bb4
begin
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using Colors
	using DataFrames
	using DataFramesMeta
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using KernelDensity
	using Latexify
	using Measurements
	using NaturalSort
	using OrderedCollections
	using Printf
	using PyCall
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
end

# ╔═╡ ebef52bc-2acf-4cf8-aca7-90cd6684c061
md"""
# Detrended Binned Light Curves

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 0158a760-1229-4089-bf90-7c7b2f1f548a
md"""
## Load data

First, let's load the relevant data needed for this notebook:
"""

# ╔═╡ 4b09c729-3395-4cee-bb69-bab59390845c
const DATA_DIR = "data/detrended/out_l_C/WASP50"

# ╔═╡ 737c135a-7412-4b87-a718-642472d4bf4b
function name(dirpath, dates_to_names)
	date_target = splitpath(split(glob(dirpath)[1], "w50_")[2])[1]
	return dates_to_names[date_target]
end

# ╔═╡ f3e9100a-ec8e-425f-9081-e457ad9b1515
dates_to_names = Dict(
	"131219_IMACS" => "Transit 1 (IMACS)",
	"150927_IMACS" => "Transit 2 (IMACS)",
	"150927_LDSS3_flat" => "Transit 2 (LDSS3)",
	"161211_IMACS" => "Transit 3 (IMACS)",
 )

# ╔═╡ 100af59b-3a24-41d0-9cda-05592bd1778f
begin
	# Load IMACS
	cubes = Dict{String, Any}()
	
	for dirpath ∈ sort(glob("$(DATA_DIR)/w50*/wavelength"))
		fpaths = sort!(glob("$(dirpath)/wbin*/PCA_1/detrended_lc.dat"), lt=natural)
		
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
		cubes[name(dirpath, dates_to_names)] = Dict(
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

# ╔═╡ b6007d1d-fb9e-4f56-a38a-febb80ea7f09
cubes |> keys

# ╔═╡ ded314ba-1ebe-4f01-bd1b-652a0258f955
@bind transit Select(cubes.keys)

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
	#const COLORS = to_colormap(:seaborn_colorblind6, 8)[[8, 6, 4, 1]]
	const COLORS = parse.(Colorant,
		[
			"#a6cee3",  # Cyan
			"#fdbf6f",  # Yellow
			"#ff7f00",  # Orange
			"#1f78b4",  # Blue
			"slategray",
			# "plum",
			# "#956cb4",  # Purple
			# "mediumaquamarine",
			# "#029e73",  # Green
		]
	)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			palette = (color=COLORS, patchcolor=COLORS,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS
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
	
	Label(fig[1:2, 0], "Relative flux + offset", rotation=π/2)
	Label(fig[end, 2:3], "Index")
	
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
		
	p |> as_svg
end

# ╔═╡ Cell order:
# ╟─ebef52bc-2acf-4cf8-aca7-90cd6684c061
# ╟─0158a760-1229-4089-bf90-7c7b2f1f548a
# ╠═4b09c729-3395-4cee-bb69-bab59390845c
# ╠═b6007d1d-fb9e-4f56-a38a-febb80ea7f09
# ╠═100af59b-3a24-41d0-9cda-05592bd1778f
# ╠═737c135a-7412-4b87-a718-642472d4bf4b
# ╠═f3e9100a-ec8e-425f-9081-e457ad9b1515
# ╟─ded314ba-1ebe-4f01-bd1b-652a0258f955
# ╠═df1a160c-22ff-4c5e-a71f-b903d8a23ef1
# ╠═bec88974-b150-4f53-9497-ddec4883ae17
# ╟─2c0406e7-96e0-4a87-a91c-02d463e32ebc
# ╠═f2da7123-cda9-47c3-aa72-4f47f4f8dfda
# ╟─1dd4968e-959a-4f6e-a0e2-9fe9b8ecdd74
# ╠═0af97a94-cb08-40e2-8011-11c8696684fa
# ╠═818783f8-7164-466e-b5a7-b75eaefe6bb4
