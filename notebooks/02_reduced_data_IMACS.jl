### A Pluto.jl notebook ###
# v0.18.0

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

# ╔═╡ b1b0690a-a1eb-11eb-1590-396d92c80c23
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	using CSV, DataFrames, DelimitedFiles, Glob, OrderedCollections
	using ImageFiltering, Statistics
	using Latexify, Printf
	using Dates
	using CondaPkg
	CondaPkg.add("numpy"); CondaPkg.resolve()
	using PythonCall
end

# ╔═╡ 9d78ea5e-e199-4509-9f74-100d1748c1c4
using JLD2

# ╔═╡ ee24f7df-c4db-4065-afe9-10be80cbcd6b
@mdx """
# Reduced data -- IMACS

In this notebook we will examine the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from IMACS.

!!! tip "Data download"
	```
	rclone sync -P ACCESS_box:WASP-50b/data/reduced_data data/reduced_data
	```
	* [Direct link](https://app.box.com/s/esq7gbpd7id98vzum1qub7twd4lba6zq)
"""

# ╔═╡ 0d766fde-8e6f-4a88-94df-49747d7c03fa
begin
	const DATA_DIR = "data/reduced_data/IMACS"
	const FIG_DIR = "figures/reduced_data"
	TableOfContents()
end

# ╔═╡ 9d180c21-e634-4a1e-8430-bdd089262f66
@mdx """
## Data extraction 🔳

The main data product from our custom pipeline for this instrument is a pickle file with the following naming scheme: `LCs_<target>_<wavelength bin scheme>.pkl`.

Each cube (`LC`) and wavelength binning scheme can be selected from the following drop-down menus, and will be used for the rest of this analysis:
"""

# ╔═╡ 28d18f7f-2e41-4771-9f27-342bbda847dd
@bind DIRPATH Select(sort(glob("$(DATA_DIR)/ut*")), default="$(DATA_DIR)/ut161211")

# ╔═╡ bd2cdf33-0c41-4948-82ab-9a28929f72b3
@bind FPATH_LC Select(glob("$(DIRPATH)/*.pkl"), default="$(DIRPATH)/LCs_w50_bins.pkl")

# ╔═╡ 3959e46c-87c9-4566-8ab1-f437323f0a9f
fname_suff = let
	suff = "IMACS_" * basename(DIRPATH)
	occursin("species", FPATH_LC) ? (suff *= "_species") : suff
end

# ╔═╡ 32b9a326-ddc8-4557-bcf5-9dcc54ed83e5
transits = merge(
	Dict("IMACS_ut$(d)" => "Transit $(i) (IMACS)"
		for (i, d) ∈ enumerate(("131219", "150927", "161211"))
	),
	Dict("IMACS_ut$(d)_species" => "Transit $(i) (IMACS) sp"
		for (i, d) ∈ enumerate(("131219", "150927", "161211"))
	),
)

# ╔═╡ e774a20f-2d58-486a-ab71-6bde678b26f8
@mdx """
## Stellar spectra ⭐

With the flux extracted for each object, we now plot the resulting stellar spectra:
"""

# ╔═╡ ac671cde-bb34-41c8-bf0d-1a52a1a0b6c6
[[1,2 ], [3, 4]] isa Array{}

# ╔═╡ 6fd88483-d005-4186-8dd2-82cea767ce90
med_std(A; dims=1) = (median(A, dims=dims), std(A, dims=dims)) .|> vec

# ╔═╡ 1f8f5bd0-20c8-4a52-9dac-4ceba18fcc06
function spec_plot!(ax, wav, A; color=:blue, norm=1.0, label="")
	μ, σ = med_std(A) ./ norm
	band!(ax, wav, μ .- σ, μ .+ σ, color=(color, 0.25))
	lines!(ax, wav, μ;
		color = color,
		cycle = Cycle(:linestyle),
		label = label,
	)
	return μ, σ
end

# ╔═╡ 818282bb-03ed-49db-9c1c-c744dad47db8
#series(LC_spectra["c28"], solid_color=:darkgrey) # Follow-up on why this forms two bands for Transit 2 (IMACS)

# ╔═╡ e3468c61-782b-4f55-a4a1-9d1883655d11
@mdx """
## White-light curves 🌅

Next, we will extract the integrated white-light curves from these spectra, divided by each comparison star:
"""

# ╔═╡ ab058d99-ce5f-4ed3-97bd-a62d2f258773
@bind window_width PlutoUI.Slider(3:2:21, default=15, show_value=true)

# ╔═╡ 7d90f304-8fc0-4b92-924c-bead3e1c0a8c
cNames_global = ("c06", "c13", "c15", "c18", "c20", "c21", "c23", "c28")

# ╔═╡ d6295509-14e1-4b24-8798-84bef7c96854
mid_transit_times = Dict(
	"Transit 1 (IMACS)" => "2013-12-19 03:22",
	"Transit 1 (IMACS) sp" => "2013-12-19 03:22",
	"Transit 2 (IMACS)" => "2015-09-27 06:37",
	"Transit 2 (IMACS) sp" => "2015-09-27 06:37",
	"Transit 3 (IMACS)" => "2016-12-12 03:03",
	"Transit 3 (IMACS) sp" => "2016-12-12 03:03",
)

# ╔═╡ 96aa3546-4ba6-4ce1-bf5c-4a00e935f702
date_fmt = dateformat"y-m-d H:M"

# ╔═╡ 06bbca64-c99f-429b-b1ad-f40f32e0deac
function compute_t_rel(t_py)
	t = pyconvert(Vector, t_py)
	t₀_utc = mid_transit_times[transits[fname_suff]]
	t₀ = DateTime(t₀_utc, date_fmt) |> datetime2julian
	return @. (t - t₀) * 24.0
end

# ╔═╡ 941cd721-07d8-4a8f-9d75-42854e6e8edb
@mdx """
!!! note
	In general, the comparison stars names (`cNames`) are not stored in alphanumeric order by default. For convenience, we ensure this sorting with `sortperm`, so that the first column corresponds to the first comparison star, the second to the second comparison star, an so on.

We next plot these light curves and identified outliers below:
"""

# ╔═╡ 4b763b58-862e-4c88-a7c9-fe0b1271c0b4
comps = Dict(
	# Transit 1
	"ut131219" => ["c15", "c18", "c21", "c23"],
	#use_comps = ["c15", "c21"]

	# Transit 2
	"ut150927" => ["c06", "c15", "c18", "c21", "c23"],
	#use_comps = ["c15", "c21"]

	# Transit 3
	"ut161211" => ["c06", "c13"],
	#use_comps = ["c06", "c13"]
)

# ╔═╡ 2df82761-b9fe-4d37-b57c-1eabb0ffa8dd
use_comps = comps[match(r"ut[0-9]{6}", fname_suff).match]

# ╔═╡ 169197fe-983d-420b-8c56-353a65b28ddc
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ 4bad8b5c-e8b9-4ceb-97f4-41b4401d4f63
@mdx """
!!! note
	We divide the target WLC by each comparison star to minimize common systematics (e.g., air mass, local refractive atmospheric effects), and to make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ ad5b07e5-75d0-4e03-a5d6-9ce4f1efd949
function filt(f_div_wlc, window_width; func=median, border="reflect")
	# Filtered light curves
	f_filt = mapslices(
		x -> mapwindow(func, x, window_width, border=border),
		f_div_wlc,
		dims = 1,
	)

	# Residuals
	Δf = f_div_wlc - f_filt

	return f_filt, abs.(Δf), Δf
end

# ╔═╡ a4517d69-76e6-462a-9449-b31d80e34a8f
# Filter specified WLCs and return superset points
function filt_idxs(f_div_wlc, window_width; ferr=0.002)
	ntimes, ncomps = size(f_div_wlc)
	f_filt, f_diff, _ = filt(f_div_wlc, window_width)
	bad_idxs = ∪(findall.(>(ferr), eachcol(f_diff))...) |> sort;
	use_idxs = deleteat!(collect(1:size(f_div_wlc, 1)), bad_idxs)
	return f_filt, use_idxs, bad_idxs
end

# ╔═╡ 0adc81ea-8678-42c2-a8b6-45fe4d26f4c4
@mdx """
### Raw flux

Just as a quick check:
"""

# ╔═╡ e6e1ea18-216a-41ae-8a1a-590793fcb669
# let
# 	fig = Figure()
# 	ax = Axis(fig[1, 1];
# 		xlabel="Index",
# 		ylabel="Relative flux",
# 		#limits=(nothing, nothing, 0.975, 1.02),
# 	)

# 	f_wlc_targ = LC_oLC ./ mean(LC_oLC, dims=1)
# 	f_wlc_comps = LC_cLC[:, sorted_cName_idxs] ./ mean(
# 		LC_cLC[:, sorted_cName_idxs], dims=1
# 	)

# 	f_wlc_comps = LC_cLC[:, sorted_cName_idxs] ./ mean(
# 		LC_cLC[:, sorted_cName_idxs], dims=1
# 	)

# 	for (i, (cName, col)) in enumerate(zip(sort(cNames), eachcol(f_wlc_comps)))
# 		scatter!(ax, col; label=cName)
# 	end

# 	scatter!(ax, f_wlc_targ ./ mean(f_wlc_targ, dims=1);
# 		linewidth = 5,
# 		color = :darkgrey,
# 		label = "WASP-50",
# 	)

# 	axislegend()

# 	fig

# end

# ╔═╡ e98dee2e-a369-448e-bfe4-8fea0f318fa8
@mdx """
## Binned light curves 🌈

We first compute `f_norm_w`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\\times`` `ncomps` ``\\times`` `nbins`, where for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve:
"""

# ╔═╡ 793c4d08-e2ee-4c9d-b7a0-11eaaddba895
@mdx """
We plot these below for each comparison star division:
"""

# ╔═╡ eeb3da97-72d5-4317-acb9-d28637a06d67
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ 06bbb3e4-9b30-43bd-941f-e357acaa80fc
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig, pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 3653ee36-35a6-4e0a-8d46-4f8389381d45
begin
	@pyexec """
	global np, pickle, load_npz, load_pickle
	import numpy as np
	import pickle

	def load_npz(fpath, allow_pickle=False):
		return np.load(fpath, allow_pickle=allow_pickle)[()]

	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	"""
	load_npz(s; allow_pickle=false) = @pyeval("load_npz")(
		s, allow_pickle=allow_pickle
	)
	load_pickle(s) = @pyeval("load_pickle")(s)
end;

# ╔═╡ dd5431a8-113c-4fa8-8fec-bf55c4b75ca4
LC = load_pickle(FPATH_LC);

# ╔═╡ 65cc9f56-1e9e-446c-82db-10dcd6334ce3
LC_spectra = pyconvert(Dict{String, Array}, LC["spectra"]);

# ╔═╡ e4e9899c-68c3-4bc4-aadd-7c5d3446e57d
save_object("/home/mango/Desktop/hmmm.jld2", LC_spectra)

# ╔═╡ bcda2043-f8c7-46bc-a5d4-b6f1f0883e9e
LC_cNames = pyconvert(Vector, LC["cNames"])

# ╔═╡ f519626c-a3e8-4390-b3af-40b7beb665ed
LC_oLC = pyconvert(Vector, LC["oLC"])

# ╔═╡ 9a9b688c-94f0-4944-a9b2-21702073e0c7
LC_cLC = pyconvert(Matrix, LC["cLC"])

# ╔═╡ 18d58341-0173-4eb1-9f01-cfa893088613
begin
	cNames = LC_cNames
	sorted_cName_idxs = sortperm(cNames)

	f_div_WLC = LC_oLC ./ LC_cLC[:, sorted_cName_idxs]
	f_div_WLC_norm = f_div_WLC ./ median(f_div_WLC, dims=1)
end

# ╔═╡ df46d106-f186-4900-9d3f-b711bc803707
@with_terminal begin
	use_comps_idxs = get_idx.(use_comps, Ref(cNames))
	_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm[:, use_comps_idxs], window_width)
	# Because python
	println(bad_idxs .- 1)
	println(use_comps_idxs .- 1)
	println(length(bad_idxs))
end

# ╔═╡ 4e4cb513-1e88-4414-aa4d-a14d934874ce
begin
	use_comps_idxs = get_idx.(use_comps, Ref(cNames))
	_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm[:, use_comps_idxs], window_width)
end;

# ╔═╡ 22b57aad-e886-4d36-bab8-baef5f3fabe6
f_div_WLC_norm

# ╔═╡ 3ca393d6-01c0-4f77-88ff-7c4f6388670e
begin
	oLCw, cLCw = pyconvert(Matrix, LC["oLCw"]), pyconvert(Array, LC["cLCw"])
	(ntimes, nbins), ncomps = size(oLCw), length(cNames)
	offs = reshape(range(0, 0.3, length=nbins), 1, :) # Arbitrary offsets for clarity
	f_norm_w = Array{Float64}(undef, ntimes, ncomps, nbins)
	for c_i ∈ 1:ncomps
		if 0.0 ∈ cLCw[:, c_i, :]
			f_w = oLCw
		else
			f_w = oLCw ./ cLCw[:, c_i, :]
		end
		f_norm_w[:, c_i, :] .= f_w ./ median(f_w, dims=1) .+ offs
	end
	baselines = ones(ntimes, nbins) .+ offs # Reference baselines
end;

# ╔═╡ a8d1c3e6-c020-495f-a443-07203b7dcd50
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
			figure_padding = 1.5,
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)

	COLORS
end

# ╔═╡ ccabf5d2-5739-4284-a972-23c02a263a5c
function plot_div_WLCS!(axs;
	t_rel, f, window_width, cNames, use_comps_idxs, ferr=0.002
)
	use_comps = cNames[use_comps_idxs]

	# Only apply filter to specified comp star divided WLCs
	f_filt, use_idxs, bad_idxs = filt_idxs(
		f[:, use_comps_idxs], window_width; ferr=ferr
	)
	
	k = 1
	c = :darkgrey
	z = 1
	for (i, cName) ∈ enumerate(cNames_global)
		if (z ≤ ncomps) && cNames_global[i] == cNames[z]
			# All points
			if cName ∈ ("c06", "c15", "c21") # LDSS3C comps
				c_text = COLORS[end]
			else
				c_text = :darkgrey
			end
			scatter!(axs[i], t_rel, f[:, z];
				color = (c, 0.3),
			)
			text!(axs[i], "$(cName)";
				position =(3, 0.98),
				align = (:right, :center),
				color = c_text,
				textsize = 24,
			)
			# Used points
			if cName ∈ use_comps
				scatter!(axs[i], t_rel[use_idxs], f[use_idxs, z];
					color = c,
				)
				lines!(axs[i], t_rel, f_filt[:, k];
					color = COLORS[end-2],
					linewidth = 2,
				)
				k += 1
			end
			z += 1 # So hacky
		else
			continue
		end

		#axislegend(axs[i])
	end
end

# ╔═╡ 13523326-a5f2-480d-9961-d23cd51192b8
let
	# Larger font for two-column
	fig = Figure(resolution=FIG_WIDE, fontsize=24)

	axs = [
		Axis(
			fig[i, j],
			limits = (-2.5, 3.5, 0.975, 1.02),
			xlabelsize = 24,
			ylabelsize = 24,
		)
		for i ∈ 1:2, j ∈ 1:4
	]
	axs = reshape(copy(fig.content), 2, 4)

	t_rel = compute_t_rel(LC["t"])
	plot_div_WLCS!(axs;
		t_rel, f=f_div_WLC_norm, window_width, cNames, use_comps_idxs
	)

	linkaxes!(axs...)
	hidexdecorations!.(axs[begin:end-1, :], grid=false)
	hideydecorations!.(axs[:, begin+1:end], grid=false)

	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2, textsize=24)
	fig[end+1, 2:end] = Label(fig, "Time from estimated mid-transit (hours)", textsize=24)

	Label(fig[0, end], transits[fname_suff];
		tellwidth = false,
		halign = :right,
		textsize = 24,
	)

	savefig(fig, "$(FIG_DIR)/div_wlcs_$(fname_suff).pdf")

	fig
end

# ╔═╡ c2c326df-1474-4b06-b183-668f0d6502aa
function plot_BLCs(datas, models, wbins, errs, comp_name; offset=0.3)
	fig = Figure(resolution=FIG_LARGE)
	median_prec = round(Int, median(errs))

	ax_left = Axis(fig[1, 1], title = "Divded BLCs")
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
			eachcol(datas),
			eachcol(models),
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

		scatter!(ax_right, baseline + resid, markersize=5, color=color)
		lines!(ax_right, baseline, linewidth=3, color=0.75*color)
		text!(ax_label, "$(wbin[1]) - $(wbin[2]) Å, $(err) ppm";
			position = (0, baseline[1]),
			textsize = 16,
			align = (:left, :center),
			offset = (-10, 2),
			color = 0.75*color,
		)
	end

	hideydecorations!.(axs[:, 2:3], grid=false)
	hidespines!(axs[end])
	hidedecorations!(axs[end])
	ylims!(ax_left, 0.95, 1.34)

	fig[1:2, 0] = Label(fig, "Relative flux + offset", rotation=π/2)
	fig[end, 2:3] = Label(fig, "Index")

	savefig(fig, "$(FIG_DIR)/div_blcs_$(fname_suff)_$(comp_name).pdf")

	fig
end

# ╔═╡ a6ec8699-475b-4a86-ab0b-d65b85de2c2d
@py begin
	import numpy as np
	import pickle
end

# ╔═╡ 6471fc66-47a5-455e-9611-c6fd9d56e9dc
wbins = pyconvert(Matrix, np.array(LC["wbins"]));

# ╔═╡ 40269026-a833-4dd8-bb22-7d26f35163e9
wbins_odd = @view wbins[begin:2:end, :] # Selects alternating bins to highlight

# ╔═╡ 589239fb-319c-40c2-af16-19025e7b28a2
begin
	# Larger font for two-column
	fig = Figure(resolution=FIG_WIDE, fontsize=24)
	ax = Axis(fig[1, 1];
		xlabel = "Wavelength (Å)",
		ylabel = "Relative flux",
		xlabelsize = 24,
		ylabelsize = 24,
		limits = (4_500, 11_000, 0.0, nothing),
	)

	wav = LC_spectra["wavelengths"]
	norm = median(LC_spectra["WASP50"])
	vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:darkgrey, 0.25))
	specs = filter(p -> p.first != "wavelengths", LC_spectra)
	
	cube_spectra = OrderedDict{String, Array}()
	cube_spectra["wav"] = wav
	cube_spectra["norm"] = [norm]
	cube_spectra["wbins"] = wbins	
	for (i, (label, f)) in enumerate(sort(specs))
		μ, σ = spec_plot!(ax, wav, f;
			color=COLORS_SERIES[i], norm, label
		)
		cube_spectra["spec_$(label)"] = [μ, σ]
	end

	axislegend(transits[fname_suff], halign=:right, gridshalign=:right)

	savefig(fig, "$(FIG_DIR)/extracted_spectra_$(fname_suff).pdf")
	save_object("/home/mango/Desktop/hmmm.jld2", cube_spectra)

	fig
end

# ╔═╡ 7962e716-8b0e-4c58-9d14-f51bbf72d419
begin
	blc_plots = OrderedDict()
	for comp_idx ∈ use_comps_idxs
		datas = f_norm_w[:, comp_idx, :]
		cName = cNames[comp_idx]
		f_med, _, f_diff = filt(datas, window_width)
		p = plot_BLCs(
			datas[use_idxs, :],
			f_med[use_idxs, :],
			wbins,
			round.(Int, reshape(std(f_diff[use_idxs, : ], dims=1), :, 1) * 1e6),
			cName,
		)
		blc_plots[cName] = p
	end
end

# ╔═╡ deb9e739-84c3-4c89-831e-1426b1ac3fbc
@bind cName Select(blc_plots.keys)

# ╔═╡ f2a2d747-0f9d-46ea-94a4-3db5b45d29c7
blc_plots[cName]

# ╔═╡ Cell order:
# ╟─ee24f7df-c4db-4065-afe9-10be80cbcd6b
# ╠═0d766fde-8e6f-4a88-94df-49747d7c03fa
# ╟─9d180c21-e634-4a1e-8430-bdd089262f66
# ╟─28d18f7f-2e41-4771-9f27-342bbda847dd
# ╟─bd2cdf33-0c41-4948-82ab-9a28929f72b3
# ╠═3959e46c-87c9-4566-8ab1-f437323f0a9f
# ╟─32b9a326-ddc8-4557-bcf5-9dcc54ed83e5
# ╠═dd5431a8-113c-4fa8-8fec-bf55c4b75ca4
# ╟─e774a20f-2d58-486a-ab71-6bde678b26f8
# ╠═589239fb-319c-40c2-af16-19025e7b28a2
# ╠═ac671cde-bb34-41c8-bf0d-1a52a1a0b6c6
# ╠═65cc9f56-1e9e-446c-82db-10dcd6334ce3
# ╠═e4e9899c-68c3-4bc4-aadd-7c5d3446e57d
# ╠═6471fc66-47a5-455e-9611-c6fd9d56e9dc
# ╠═40269026-a833-4dd8-bb22-7d26f35163e9
# ╠═1f8f5bd0-20c8-4a52-9dac-4ceba18fcc06
# ╠═6fd88483-d005-4186-8dd2-82cea767ce90
# ╠═818282bb-03ed-49db-9c1c-c744dad47db8
# ╟─e3468c61-782b-4f55-a4a1-9d1883655d11
# ╟─ab058d99-ce5f-4ed3-97bd-a62d2f258773
# ╠═13523326-a5f2-480d-9961-d23cd51192b8
# ╠═7d90f304-8fc0-4b92-924c-bead3e1c0a8c
# ╠═ccabf5d2-5739-4284-a972-23c02a263a5c
# ╠═d6295509-14e1-4b24-8798-84bef7c96854
# ╠═96aa3546-4ba6-4ce1-bf5c-4a00e935f702
# ╠═06bbca64-c99f-429b-b1ad-f40f32e0deac
# ╠═bcda2043-f8c7-46bc-a5d4-b6f1f0883e9e
# ╠═f519626c-a3e8-4390-b3af-40b7beb665ed
# ╠═9a9b688c-94f0-4944-a9b2-21702073e0c7
# ╠═18d58341-0173-4eb1-9f01-cfa893088613
# ╟─941cd721-07d8-4a8f-9d75-42854e6e8edb
# ╠═4b763b58-862e-4c88-a7c9-fe0b1271c0b4
# ╠═2df82761-b9fe-4d37-b57c-1eabb0ffa8dd
# ╠═df46d106-f186-4900-9d3f-b711bc803707
# ╠═169197fe-983d-420b-8c56-353a65b28ddc
# ╟─4bad8b5c-e8b9-4ceb-97f4-41b4401d4f63
# ╠═4e4cb513-1e88-4414-aa4d-a14d934874ce
# ╠═22b57aad-e886-4d36-bab8-baef5f3fabe6
# ╠═a4517d69-76e6-462a-9449-b31d80e34a8f
# ╠═ad5b07e5-75d0-4e03-a5d6-9ce4f1efd949
# ╟─0adc81ea-8678-42c2-a8b6-45fe4d26f4c4
# ╠═e6e1ea18-216a-41ae-8a1a-590793fcb669
# ╟─e98dee2e-a369-448e-bfe4-8fea0f318fa8
# ╟─deb9e739-84c3-4c89-831e-1426b1ac3fbc
# ╟─f2a2d747-0f9d-46ea-94a4-3db5b45d29c7
# ╟─7962e716-8b0e-4c58-9d14-f51bbf72d419
# ╠═3ca393d6-01c0-4f77-88ff-7c4f6388670e
# ╟─793c4d08-e2ee-4c9d-b7a0-11eaaddba895
# ╠═c2c326df-1474-4b06-b183-668f0d6502aa
# ╟─eeb3da97-72d5-4317-acb9-d28637a06d67
# ╟─06bbb3e4-9b30-43bd-941f-e357acaa80fc
# ╠═3653ee36-35a6-4e0a-8d46-4f8389381d45
# ╠═a8d1c3e6-c020-495f-a443-07203b7dcd50
# ╠═a6ec8699-475b-4a86-ab0b-d65b85de2c2d
# ╠═b1b0690a-a1eb-11eb-1590-396d92c80c23
# ╠═9d78ea5e-e199-4509-9f74-100d1748c1c4
