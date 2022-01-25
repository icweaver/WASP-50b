### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 818783f8-7164-466e-b5a7-b75eaefe6bb4
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	using AlgebraOfGraphics
	using CairoMakie
	import CairoMakie.Makie.KernelDensity: kde
	using CSV, DataFrames
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using Measurements, Statistics
	using OrderedCollections
	using Printf
	using NaturalSort
end

# ╔═╡ ebef52bc-2acf-4cf8-aca7-90cd6684c061
md"""
# Detrended binned light curves

In this notebook we will plot the detrended binned light curves for all nights.

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 6008853d-1a74-4004-8aa8-7a70d8297045
md"""
!!! note "Data download"
	```shell
	rsync -azRP $H:/pool/sao_access/iweaver/GPTransmissionSpectra/./"out_*" data/detrended/ --exclude={"*george*","*mnest*"}
	```

	```shell
	rclone sync -P ACCESS_box:WASP-50b/data/detrended data/detrended
	```
"""

# ╔═╡ 0158a760-1229-4089-bf90-7c7b2f1f548a
md"""
## Load data

First, let's load the relevant data needed for this notebook:
"""

# ╔═╡ 4b09c729-3395-4cee-bb69-bab59390845c
@bind DATA_DIR Select(glob("data/detrended/out_*/WASP50"))

# ╔═╡ b63189de-7f78-46d3-a119-3064e275dbe4
const FIG_PATH = "figures/detrended"

# ╔═╡ 737c135a-7412-4b87-a718-642472d4bf4b
function name(dirpath, dates_to_names)
	date_target = splitpath(split(glob(dirpath)[1], "w50_")[2])[1]
	return dates_to_names[date_target]
end

# ╔═╡ f3e9100a-ec8e-425f-9081-e457ad9b1515
dates_to_names = Dict(
	"131219_IMACS" => "Transit_1_IMACS",
	"131219_sp_IMACS" => "Transit_1_IMACS_sp",
	"150927_IMACS" => "Transit_2_IMACS",
	"150927_sp_IMACS" => "Transit_2_IMACS_sp",
	"150927_LDSS3" => "Transit_2_LDSS3",
	"150927_sp_LDSS3" => "Transit_2_LDSS3_sp",
	"161211_IMACS" => "Transit_3_IMACS",
	"161211_sp_IMACS" => "Transit_3_IMACS_sp"
 )

# ╔═╡ 100af59b-3a24-41d0-9cda-05592bd1778f
begin
	# Load IMACS
	cubes = Dict{String, Any}()
	
	for dirpath ∈ sort(glob("$(DATA_DIR)/w50*/wavelength"))
		fpaths = sort!(glob("$(dirpath)/wbin*/PCA_1/detrended_lc.dat"), lt=natural)
		dirpath_WLC = "$(dirname(dirpath))/white-light"

		# TODO, track this down
		deleteat!(fpaths, findfirst(s -> occursin("wbin3", s), fpaths))
		
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
			df = CSV.read(fpath, DataFrame;
				header=["Time", "DetFlux", "DetFluxErr", "Model"],
				comment = "#",
				select=[:DetFlux, :Model],
			)
			#@show fpath nrow(df)
			lc .= df.DetFlux
			model .= df.Model
		end
		
		# Save
		transit = name(dirpath, dates_to_names)
		cubes[transit] = Dict(
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

# ╔═╡ a5acf744-dfe7-4088-885b-9142af8f0d8f
@bind transit Select(cubes.keys)

# ╔═╡ 1dd4968e-959a-4f6e-a0e2-9fe9b8ecdd74
md"""
## Notebook setup
"""

# ╔═╡ c59e2697-d2a3-4bdb-ba64-059246697c1c
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 0af97a94-cb08-40e2-8011-11c8696684fa
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (800, 600)
	const FIG_LARGE = (1_200, 1_000)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#a6cee3",  # Cyan
			"#fdbf6f",  # Yellow
			"#ff7f00",  # Orange
			"#1f78b4",  # Blue
		]
	)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (
				xlabelsize = 18,
				ylabelsize = 18,
				topspinevisible = true,
				rightspinevisible = true,
				topspinecolor = :darkgrey,
				rightspinecolor = :darkgrey
			),
			Label = (textsize=18,  padding=(0, 10, 0, 0)),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)
	
	COLORS
end

# ╔═╡ bec88974-b150-4f53-9497-ddec4883ae17
function plot_BLCs(transit, datas, models, wbins, errs; offset=0.3)
	fig = Figure(resolution=FIG_TALL)
	median_prec = round(Int, median(errs))
	ax_left = Axis(fig[1, 1], title = "Detrended BLCs")
	ax_right = Axis(fig[1, 2], title = "Residuals")
	ax_label = Axis(fig[1, 3], title = "Median precision: $(median_prec) ppm")
	axs = reshape(copy(fig.content), 1, 3)
	#ylims!(ax_label, (1.05, 1.15))
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
	ylims!(ax_left, 0.95, 1.34)
	
	Label(fig[1:2, 0], "Relative flux + offset", rotation=π/2)
	Label(fig[end, 2:3], "Index")

	f_suff = basename(dirname(DATA_DIR))
	savefig(fig, "$(FIG_PATH)/detrended_blcs_$(transit)_$(f_suff).png")
	
	fig
end

# ╔═╡ df1a160c-22ff-4c5e-a71f-b903d8a23ef1
begin
	blc_plots = Dict()
	for transit in cubes.keys
		tspec = cubes[transit]["tspec"]
		errs = maximum([tspec.Depthup_ppm_ tspec.DepthDown_ppm_], dims=2)
		
		p = plot_BLCs(
			transit,
			cubes[transit]["fluxes"],
			cubes[transit]["models"],
			cubes[transit]["wbins"],
			round.(Int, errs),
		)
		blc_plots[transit] = p
	end
end

# ╔═╡ cb2a3117-03be-42f4-adf6-c23a42252ddf
blc_plots[transit]

# ╔═╡ Cell order:
# ╟─ebef52bc-2acf-4cf8-aca7-90cd6684c061
# ╟─6008853d-1a74-4004-8aa8-7a70d8297045
# ╟─0158a760-1229-4089-bf90-7c7b2f1f548a
# ╟─4b09c729-3395-4cee-bb69-bab59390845c
# ╟─b63189de-7f78-46d3-a119-3064e275dbe4
# ╠═100af59b-3a24-41d0-9cda-05592bd1778f
# ╠═737c135a-7412-4b87-a718-642472d4bf4b
# ╠═f3e9100a-ec8e-425f-9081-e457ad9b1515
# ╠═df1a160c-22ff-4c5e-a71f-b903d8a23ef1
# ╟─a5acf744-dfe7-4088-885b-9142af8f0d8f
# ╟─cb2a3117-03be-42f4-adf6-c23a42252ddf
# ╠═bec88974-b150-4f53-9497-ddec4883ae17
# ╟─1dd4968e-959a-4f6e-a0e2-9fe9b8ecdd74
# ╟─c59e2697-d2a3-4bdb-ba64-059246697c1c
# ╠═0af97a94-cb08-40e2-8011-11c8696684fa
# ╠═818783f8-7164-466e-b5a7-b75eaefe6bb4
