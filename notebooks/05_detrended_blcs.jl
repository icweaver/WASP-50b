### A Pluto.jl notebook ###
# v0.18.1

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

# ╔═╡ d5ffca19-ab35-41f8-81f0-a2972de7f758
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	import CairoMakie.Makie.KernelDensity: kde
	using CSV, DataFrames, DataFramesMeta, DelimitedFiles, Glob, OrderedCollections
	using ImageFiltering, Measurements, Statistics
	using Latexify, Printf
	using Dates, NaturalSort
end

# ╔═╡ 6c6741b8-eeb1-4c1e-8d22-40d08df00ced
begin
	const BASE_DIR = "data/detrended"
	const FIG_DIR = "figures/detrended"
	TableOfContents()
end

# ╔═╡ ebef52bc-2acf-4cf8-aca7-90cd6684c061
@mdx """
# Detrended binned light curves

In this notebook we will plot the detrended binned light curves for all nights.

!!! tip "Data download"
	```
	rclone sync -P ACCESS_box:WASP-50b/$(BASE_DIR) $(BASE_DIR)
	```

	* [Direct link](https://app.box.com/s/wr8tpof238cq8oj71ulaf69z9q0k7f9w)
"""

# ╔═╡ 0158a760-1229-4089-bf90-7c7b2f1f548a
@mdx """
## Load data ⬇️

First, let's load the relevant data needed for this notebook:
"""

# ╔═╡ 4b09c729-3395-4cee-bb69-bab59390845c
@bind DATA_DIR Select(glob("$(BASE_DIR)/out_*/WASP50"))

# ╔═╡ 141a652d-0d43-4ebc-9a8c-ac43f31e7831
fname_suff = (basename ∘ dirname)(DATA_DIR)

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
	"150927_LDSS3_flat" => "Transit_2_LDSS3",
	"150927_sp_LDSS3_flat" => "Transit_2_LDSS3_sp",
	"161211_IMACS" => "Transit_3_IMACS",
	"161211_sp_IMACS" => "Transit_3_IMACS_sp"
 )

# ╔═╡ 100af59b-3a24-41d0-9cda-05592bd1778f
begin
	# Load IMACS
	cubes = Dict{String, Any}()

	for dirpath ∈ sort(glob("$(DATA_DIR)/w50*/wavelength"))
		fpaths = sort!(glob("$(dirpath)/wbin*/PCA_1/detrended_lc.dat"), lt=natural)
		#@show fpaths
		dirpath_WLC = "$(dirname(dirpath))/white-light"

		# TODO, track this down
		#deleteat!(fpaths, findfirst(s -> occursin("wbin11", s), fpaths))

		# WLC BMA t₀ and P
		t₀, P = let
			df_BMA_WLC = CSV.read("$(dirpath_WLC)/results.dat", DataFrame;
				comment = "#",
				normalizenames = true,
				stripwhitespace = true,
			)
			(df_BMA_WLC[findfirst(==("t0"), df_BMA_WLC.Variable), :].Value,
			df_BMA_WLC[findfirst(==("P"), df_BMA_WLC.Variable), :].Value)
		end

		# Pre-allocate matrices to hold detrended data and models in each bin
		N_time = size(readdlm("$(dirpath_WLC)/comps.dat", ','), 1)
		N_bins = length(fpaths)
		det_BLC_fluxes = Matrix{Float64}(undef, N_time, N_bins)
		det_BLC_models = copy(det_BLC_fluxes)
		# Populate each matrix
		times = []
		for (fpath, lc, model) in zip(
			fpaths, eachcol(det_BLC_fluxes), eachcol(det_BLC_models)
		)
			df = CSV.read(fpath, DataFrame;
				header=["Time", "DetFlux", "DetFluxErr", "Model"],
				comment = "#",
			)
			#@show fpath nrow(df)
			lc .= df.DetFlux
			model .= df.Model
			push!(times, df.Time)
		end

		# Save
		transit = name(dirpath, dates_to_names)
		fpath = "$(dirname(dirpath_WLC))/transpec.csv"
		cubes[transit] = Dict(
			"fluxes" => det_BLC_fluxes,
			"models" => det_BLC_models,
			"t" => times[begin],
			"t₀" => t₀,
			"P" => P,
			"tspec" => CSV.read(
				fpath, DataFrame;
				normalizenames = true,
			)
		)
		# Add wbins for LDSS3C
		df = cubes[transit]["tspec"]
		if occursin("_sp_LDSS3", fpath)
			df.Wav_d, df.Wav_u = eachcol(readdlm("data/detrended/wbins/w50_bins_species.dat"; comments=true))
		elseif occursin("_LDSS3", fpath)
			df.Wav_d, df.Wav_u = eachcol(readdlm("data/detrended/wbins/w50_bins_LDSS3.dat"; comments=true))
		end
	end

	cubes = sort(cubes)
end

# ╔═╡ 03f34fcd-61e0-4996-8e9e-19c3e0bb1e88
cubes["Transit_3_IMACS"]["t"][[begin, end]] .|> julian2datetime .|> (x -> round(x, Minute))

# ╔═╡ efb8ed46-1607-4c13-b1ba-e4ca37e59b98
@mdx """
## Plot 📊
"""

# ╔═╡ a5acf744-dfe7-4088-885b-9142af8f0d8f
@bind transit Select(cubes.keys)

# ╔═╡ ff8fb667-fcaf-4c2d-b5c9-16252bcfeade
# Computes orbital phase
function compute_ϕ(t, t₀, P)
	phase = ((t - t₀) / P) % 1.0
    phase ≥ 0.5 && (phase -= 1.0)
	return phase
end

# ╔═╡ 1dd4968e-959a-4f6e-a0e2-9fe9b8ecdd74
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ c59e2697-d2a3-4bdb-ba64-059246697c1c
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig, pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 0af97a94-cb08-40e2-8011-11c8696684fa
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = 72 .* (6, 8)
	const FIG_WIDE = 72 .* (12, 6)
	const FIG_LARGE = 72 .* (12, 12)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#66C2A5",  # Green
			"#FDBF6F",  # Yellow
			"#FF7F00",  # Orange
			"#1F78B4",  # Blue
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
				rightspinecolor = :darkgrey,
			),
			Label = (
				textsize = 18,
				font = AlgebraOfGraphics.firasans("Medium"),
			),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Regular"), textsize=18),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			figure_padding = (0, 1.5, 0, 0),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)

	COLORS
end

# ╔═╡ bec88974-b150-4f53-9497-ddec4883ae17
function plot_BLCs(transit, phase, datas, models, errs, wbins; offset=0.3)
	fig = Figure(resolution=FIG_LARGE)
	median_prec = round(Int, median(errs))
	ax_left = Axis(fig[1, 1], title = "Detrended BLCs")
	ax_right = Axis(fig[1, 2], title = "Residuals")
	ax_label = Axis(fig[1, 3], title = "Median precision: $(median_prec) ppm")
	axs = reshape(copy(fig.content), 1, 3)
	
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
		scatter!(ax_left, phase, data, strokewidth=0, markersize=5, color=color)
		lines!(ax_left, phase, model, linewidth=3, color=0.75*color)
		
		scatter!(ax_right, phase, baseline + resid, markersize=5, color=color)
		lines!(ax_right, phase, baseline, linewidth=3, color=0.75*color)
		
		text!(ax_label, "$(wbin[1]) - $(wbin[2]) Å, $(err) ppm";
			position = (-0.03, baseline[1]),
			textsize = 16,
			align = (:left, :center),
			#offset = (-10, 2),
			color = 0.75*color,
		)
	end
	
	hideydecorations!.(axs[:, 2:3], grid=false)
	hidespines!(axs[end])
	hidedecorations!(axs[end])
	xlims!(ax_left, -0.04, 0.065)
	ylims!(ax_left, 0.95, 1.34)
	linkaxes!(axs...)
	
	Label(fig[1:2, 0], "Relative flux + offset", rotation=π/2)
	Label(fig[end, 2:3], "Phase")

	#savefig(fig, "$(FIG_DIR)/detrended_blcs_$(transit)_$(fname_suff).pdf")
	savefig(fig, "/home/mango/Desktop/detrended_blcs_$(transit)_$(fname_suff).png")
	
	fig
end

# ╔═╡ df1a160c-22ff-4c5e-a71f-b903d8a23ef1
begin
	blc_plots = Dict()
	for (transit, cube) ∈ cubes
		tspec = cube["tspec"]
		wbins = tspec[:, [:Wav_d, :Wav_u]]
		errs = mean([tspec.Depthup_ppm_ tspec.DepthDown_ppm_], dims=2)
		t, t₀, P = get.(Ref(cube), ["t", "t₀", "P"], nothing)
		phase = compute_ϕ.(t, t₀, P)
		
		p = plot_BLCs(
			transit,
			phase,
			cube["fluxes"],
			cube["models"],
			round.(Int, errs),
			wbins,
		)
		blc_plots[transit] = p
	end
end

# ╔═╡ cb2a3117-03be-42f4-adf6-c23a42252ddf
blc_plots[transit]

# ╔═╡ Cell order:
# ╟─ebef52bc-2acf-4cf8-aca7-90cd6684c061
# ╠═6c6741b8-eeb1-4c1e-8d22-40d08df00ced
# ╟─0158a760-1229-4089-bf90-7c7b2f1f548a
# ╟─4b09c729-3395-4cee-bb69-bab59390845c
# ╠═100af59b-3a24-41d0-9cda-05592bd1778f
# ╠═141a652d-0d43-4ebc-9a8c-ac43f31e7831
# ╠═737c135a-7412-4b87-a718-642472d4bf4b
# ╠═f3e9100a-ec8e-425f-9081-e457ad9b1515
# ╠═03f34fcd-61e0-4996-8e9e-19c3e0bb1e88
# ╟─efb8ed46-1607-4c13-b1ba-e4ca37e59b98
# ╠═df1a160c-22ff-4c5e-a71f-b903d8a23ef1
# ╟─a5acf744-dfe7-4088-885b-9142af8f0d8f
# ╟─cb2a3117-03be-42f4-adf6-c23a42252ddf
# ╠═bec88974-b150-4f53-9497-ddec4883ae17
# ╠═ff8fb667-fcaf-4c2d-b5c9-16252bcfeade
# ╟─1dd4968e-959a-4f6e-a0e2-9fe9b8ecdd74
# ╟─c59e2697-d2a3-4bdb-ba64-059246697c1c
# ╠═0af97a94-cb08-40e2-8011-11c8696684fa
# ╠═d5ffca19-ab35-41f8-81f0-a2972de7f758
