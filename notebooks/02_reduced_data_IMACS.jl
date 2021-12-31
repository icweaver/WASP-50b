### A Pluto.jl notebook ###
# v0.17.3

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
	using AlgebraOfGraphics
	using CairoMakie
	using DataFrames
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using Statistics
end

# ╔═╡ d07e895f-c43b-46c4-be47-5af44bfda47a
begin
	ENV["PYTHON"] = ""
	Pkg.build("PyCall")
	using PyCall
	const Conda = PyCall.Conda
	Conda.add("lightkurve", :WASP50b)
end

# ╔═╡ ee24f7df-c4db-4065-afe9-10be80cbcd6b
md"""
# Reduced data -- IMACS

In this notebook we will examine the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from IMACS.

$(TableOfContents(depth=4))
"""

# ╔═╡ 9d180c21-e634-4a1e-8430-bdd089262f66
md"""
## Data extraction 🔳

The main data product from our custom pipeline for this instrument is a pickle file with the following naming scheme: `LCs_<target>_<wavelength bin scheme>.pkl`. 

Each cube (`LC`) can be selected from the following drop-down menus, and will be used for the rest of this analysis:
"""

# ╔═╡ 28d18f7f-2e41-4771-9f27-342bbda847dd
@bind DIRPATH Select(sort(glob("data/reduced_data/IMACS/ut*")))

# ╔═╡ bd2cdf33-0c41-4948-82ab-9a28929f72b3
@bind FPATH_LC Select(glob("$(DIRPATH)/*.pkl"))

# ╔═╡ 3959e46c-87c9-4566-8ab1-f437323f0a9f
fname_suff = let
	suff = basename(DIRPATH)
	occursin("species", FPATH_LC) ? (suff *= "_species") : suff
end

# ╔═╡ 6471fc66-47a5-455e-9611-c6fd9d56e9dc
wbins = let
	if occursin("species", FPATH_LC)
		wbins_name = "w50_bins_species.dat"
	elseif occursin("ut131219", FPATH_LC)
		wbins_name = "w50_bins_ut131219.dat"
	else
		wbins_name = "w50_bins.dat"
	end
	readdlm("$(DIRPATH)/../../wbins/$(wbins_name)", comments=true)
end

# ╔═╡ 7f30d3c2-d8a4-4c0d-bbca-b1d5b1e4c20b
with_terminal() do
	wl, wu = wbins[:, 1], wbins[:, 2]
	Δw = @. wu - wl
	wc = @. (wl + wu) / 2
	df = DataFrame(
		"Central wavelength" => wc,
		"Lower wav." => wl,
		"Upper wav." => wu,
		L"$\Delta$ wav." => Δw,
	)
	latextabular(df, latex=false) |> println
end

# ╔═╡ 66052b03-35a0-4877-abef-f525766fd332
md"""
!!! tip
	A description of each field can be found on our repo page here (PUBLIC LINK?)
"""

# ╔═╡ 7bfc971c-8737-49ad-adec-ac57d176f10e
md"""
We can extract the comparison star flux in a similar way by stacking the ``N \times W`` matrix for each star:
"""

# ╔═╡ 1c3e8cb3-2eff-47c2-8c17-01d0599556b8
md"""
### Helper functions
"""

# ╔═╡ 3653ee36-35a6-4e0a-8d46-4f8389381d45
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

# ╔═╡ dd5431a8-113c-4fa8-8fec-bf55c4b75ca4
LC = load_pickle(FPATH_LC);

# ╔═╡ e774a20f-2d58-486a-ab71-6bde678b26f8
md"""
## Stellar spectra ⭐

With the flux extracted for each object, we now turn to analyzing the resulting stellar spectra:
"""

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
end

# ╔═╡ e3468c61-782b-4f55-a4a1-9d1883655d11
md"""
## White-light curves 🌅

Next, we will extract the integrated white-light curves from these spectra, divided by each comparison star:
"""

# ╔═╡ bcda2043-f8c7-46bc-a5d4-b6f1f0883e9e
LC_cNames = LC["cNames"]

# ╔═╡ f519626c-a3e8-4390-b3af-40b7beb665ed
LC_oLC = LC["oLC"]

# ╔═╡ 9a9b688c-94f0-4944-a9b2-21702073e0c7
LC_cLC = LC["cLC"]

# ╔═╡ 18d58341-0173-4eb1-9f01-cfa893088613
begin
	comp_names = LC_cNames
	sorted_cName_idxs = sortperm(comp_names)
	
	f_div_WLC = LC_oLC ./ LC_cLC[:, sorted_cName_idxs]
	f_div_WLC_norm = f_div_WLC ./ median(f_div_WLC, dims=1)
end

# ╔═╡ 941cd721-07d8-4a8f-9d75-42854e6e8edb
md"""
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
	"ut150927" => ["c15", "c18", "c21", "c23"],
	#use_comps = ["c15", "c21"]

	# Transit 3
	"ut161211" => ["c06", "c13"],
	#use_comps = ["c06", "c13"]
)

# ╔═╡ 2df82761-b9fe-4d37-b57c-1eabb0ffa8dd
use_comps = comps[split(fname_suff, '_')[1]]

# ╔═╡ ab058d99-ce5f-4ed3-97bd-a62d2f258773
@bind window_width PlutoUI.Slider(3:2:21, default=15, show_value=true)

# ╔═╡ 169197fe-983d-420b-8c56-353a65b28ddc
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ 4bad8b5c-e8b9-4ceb-97f4-41b4401d4f63
md"""
!!! note
	We divide the target WLC by each comparison star to minimize common systematics (e.g., air mass, local refractive atmospheric effects), and to make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ 22b57aad-e886-4d36-bab8-baef5f3fabe6
f_div_WLC_norm

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

# ╔═╡ df46d106-f186-4900-9d3f-b711bc803707
@with_terminal begin
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names))
	_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm[:, use_comps_idxs], window_width)
	# Because python
	println(bad_idxs .- 1)
	println(use_comps_idxs .- 1)
	println(length(bad_idxs))
end

# ╔═╡ 4e4cb513-1e88-4414-aa4d-a14d934874ce
begin
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names))
	_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm[:, use_comps_idxs], window_width)
end;

# ╔═╡ 0adc81ea-8678-42c2-a8b6-45fe4d26f4c4
md"""
### Raw flux
"""

# ╔═╡ e6e1ea18-216a-41ae-8a1a-590793fcb669
let
	fig = Figure()
	ax = Axis(fig[1, 1];
		xlabel="Index",
		ylabel="Relative flux",
		limits=(nothing, nothing, 0.975, 1.02),
	)
	
	f_wlc_targ = LC_oLC ./ mean(LC_oLC, dims=1)
	f_wlc_comps = LC_cLC[:, sorted_cName_idxs] ./ mean(
		LC_cLC[:, sorted_cName_idxs], dims=1
	)
	
	f_wlc_comps = LC_cLC[:, sorted_cName_idxs] ./ mean(
		LC_cLC[:, sorted_cName_idxs], dims=1
	)
	
	for (i, (cName, col)) in enumerate(zip(sort(comp_names), eachcol(f_wlc_comps)))
		scatter!(ax, col; label=cName)
	end
	
	scatter!(ax, f_wlc_targ ./ mean(f_wlc_targ, dims=1);
		linewidth = 5,
		color = :darkgrey,
		label = "WASP-50",
	)
	
	axislegend()
	
	fig
		
end

# ╔═╡ e98dee2e-a369-448e-bfe4-8fea0f318fa8
md"""
## Binned light curves 🌈

We first compute `f_norm_w`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\times`` `ncomps` ``\times`` `nbins`, where for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve:
"""

# ╔═╡ 793c4d08-e2ee-4c9d-b7a0-11eaaddba895
md"""
We plot these below for each comparison star division. Move the slider to view the plot for the corresponding comparison star:

comp star $(@bind comp_idx PlutoUI.Slider(use_comps_idxs))

!!! note "Future"
	Ability to interact with sliders completely in the browser coming soon!
"""

# ╔═╡ 3ca393d6-01c0-4f77-88ff-7c4f6388670e
begin
	oLCw, cLCw = LC["oLCw"], LC["cLCw"]
	(ntimes, nbins), ncomps = size(oLCw), length(comp_names)
	offs = reshape(range(0, 0.3, length=nbins), 1, :) # Arbitrary offsets for clarity
	f_norm_w = Array{Float64}(undef, ntimes, ncomps, nbins)
	for c_i in 1:ncomps
		f_w = oLCw ./ cLCw[:, c_i, :]
		f_norm_w[:, c_i, :] .= f_w ./ median(f_w, dims=1) .+ offs
	end
	baselines = ones(size(f_norm_w[:, comp_idx, :])) .+ offs # Reference baselines
end;

# ╔═╡ 2768623f-7904-4674-a2ee-ad809cdd508b
# Median filtered curves for visualization purposes
f_med, _, f_diff = filt(f_norm_w[:, comp_idx, :], window_width)

# ╔═╡ eeb3da97-72d5-4317-acb9-d28637a06d67
md"""
## Notebook setup
"""

# ╔═╡ a8d1c3e6-c020-495f-a443-07203b7dcd50
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (800, 600)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#a6cee3",  # Cyan
			"#fdbf6f",  # Yellow
			"#ff7f00",  # Orange
			"#1f78b4",  # Blue
			# "plum",
			# "#956cb4",  # Purple
			# "mediumaquamarine",
			# "#029e73",  # Green,
		]
	)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,  padding=(0, 10, 0, 0)),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS
end

# ╔═╡ 589239fb-319c-40c2-af16-19025e7b28a2
let
	fig = Figure(resolution=FIG_WIDE)
	ax = Axis(fig[1, 1], xlabel="Wavelength (Å)", ylabel="Relative flux")
	
	LC_spectra = LC["spectra"]
	wav = LC_spectra["wavelengths"]
	f_norm = median(LC_spectra["WASP50"])
	
	i = 1
	for (name, f) in sort(LC_spectra)
		if name != "wavelengths"
			spec_plot!(ax, wav, f, color=COLORS_SERIES[i], norm=f_norm, label=name)
			i += 1
		end
	end
	
	vlines!.(ax, wbins, linewidth=1.0, color=:lightgrey)
	
	axislegend()
	
	xlims!(ax, 4_500, 11_000)
	ylims!(ax, 0, 2.6)
	
	path = "../../ACCESS_WASP-50b/figures/reduced_data"
	mkpath(path)
	save("$(path)/extracted_spectra_IMACS_$(fname_suff).png", fig)
	
	fig
end

# ╔═╡ ccabf5d2-5739-4284-a972-23c02a263a5c
function plot_div_WLCS!(
	axs, f_div_wlc, window_width, cNames, use_comps_idxs; ferr=0.002
)
	use_comps = cNames[use_comps_idxs]
	
	# Only apply filter to specified comp star divided WLCs
	f_filt, use_idxs, bad_idxs = filt_idxs(
		f_div_wlc[:, use_comps_idxs], window_width; ferr=ferr
	)
	idxs = 1:size(f_div_wlc, 1)
	k = 1
	c = :darkgrey
	for (i, cName) ∈ enumerate(cNames)		
		# All points
		if cName ∈ ("c06", "c15", "c21") # LDSS3 comps
			c_text = COLORS[end]
		else
			c_text = :darkgrey
		end
		scatter!(axs[i], idxs, f_div_wlc[:, i];
			color = (c, 0.3),
		)
		text!(axs[i], "$(cName)";
			position =(300, 0.975),
			align = (:right, :center),
			color = c_text,
		)
		
		# Used points
		if cName ∈ use_comps
			scatter!(axs[i], idxs[use_idxs], f_div_wlc[use_idxs, i];
				color = c,
			)
			lines!(axs[i], idxs, f_filt[:, k];
				color = COLORS[end-2],
				linewidth = 2,
			)
			k += 1
		end
		
		#axislegend(axs[i])
	end
end

# ╔═╡ 13523326-a5f2-480d-9961-d23cd51192b8
let
	fig = Figure(resolution=FIG_WIDE)
	
	axs = [Axis(fig[i, j]) for i ∈ 1:2, j ∈ 1:4]
	axs = reshape(copy(fig.content), 2, 4)
	
	plot_div_WLCS!(
		axs, f_div_WLC_norm, window_width, comp_names, use_comps_idxs
	)
	
	linkaxes!(axs...)
	hidexdecorations!.(axs[begin:end-1, :], grid=false)
	hideydecorations!.(axs[:, begin+1:end], grid=false)
	ylims!(axs[end], 0.97, 1.02)
	
	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2)
	fig[end+1, 2:end] = Label(fig, "Index")
	
	path = "../../ACCESS_WASP-50b/figures/reduced_data"
	mkpath(path)
	save("$(path)/div_wlcs_IMACS_$(fname_suff).png", fig)
	
	fig
end

# ╔═╡ 684c026a-b5d0-4694-8d29-a44b7cb0fd6c
let
	fig = Figure(resolution=FIG_TALL)
	
	comp_name = comp_names[comp_idx]
	ax_left = Axis(fig[1, 1], title = "target / $(comp_name)")
	ax_right = Axis(fig[1, 2], title = "residuals")
	ax_label = Axis(fig[1, 3])
	axs = reshape(copy(fig.content), 1, 3)
	linkaxes!(axs...)
		
	colors = to_colormap(:Spectral_4, nbins) |> reverse
	for (f, f_med, resid, b, c, w) in zip(
			eachcol(f_norm_w[:, comp_idx, :]),
			eachcol(f_med),
			eachcol(f_diff),
			eachcol(baselines),
			colors,
			eachrow(wbins),
		)	
		scatter!(ax_left, f[use_idxs], markersize=5, color=c)
		lines!(ax_left, f_med[use_idxs], linewidth=3, color=0.75*c)
		
		scatter!(ax_right, (b + resid)[use_idxs], markersize=5, color=c)
		lines!(ax_right, b[use_idxs], linewidth=3, color=0.75*c)
		text!(ax_label, "$(w[1]) - $(w[2]) Å";
			position = Point2f0(0, b[1]),
			textsize = 16,
			align = (:left, :center),
			offset = Point2f0(0, 2),
			color = 0.75*c,
		)
	end
	
	hideydecorations!.(axs[:, 2:3], grid=false)
	hidespines!(axs[end])
	hidedecorations!(axs[end])
	ylims!(ax_left, 0.95, 1.34)
	
	fig[1:2, 0] = Label(fig, "Relative flux + offset", rotation=π/2)
	fig[end, 2:3] = Label(fig, "Index")
	
	path = "../../ACCESS_WASP-50b/figures/reduced_data"
	mkpath(path)
	save("$(path)/div_blcs_IMACS_$(fname_suff)_$(comp_name).png", fig)
	
	fig
end

# ╔═╡ e8478e36-10fc-4e95-bf09-e217cad0cb15
@with_terminal Conda.list(:WASP50b)

# ╔═╡ 03af71ac-673b-459b-a931-a600b13d7ee6
html"""
<style>
#launch_binder {
	display: none;
}
body.disable_ui main {
		max-width : 95%;
	}
@media screen and (min-width: 1081px) {
	body.disable_ui main {
		margin-left : 10px;
		max-width : 72%;
		align-self: flex-start;
	}
}
</style>
"""

# ╔═╡ Cell order:
# ╟─ee24f7df-c4db-4065-afe9-10be80cbcd6b
# ╟─9d180c21-e634-4a1e-8430-bdd089262f66
# ╟─28d18f7f-2e41-4771-9f27-342bbda847dd
# ╟─bd2cdf33-0c41-4948-82ab-9a28929f72b3
# ╠═3959e46c-87c9-4566-8ab1-f437323f0a9f
# ╠═dd5431a8-113c-4fa8-8fec-bf55c4b75ca4
# ╟─6471fc66-47a5-455e-9611-c6fd9d56e9dc
# ╟─7f30d3c2-d8a4-4c0d-bbca-b1d5b1e4c20b
# ╟─66052b03-35a0-4877-abef-f525766fd332
# ╟─7bfc971c-8737-49ad-adec-ac57d176f10e
# ╟─1c3e8cb3-2eff-47c2-8c17-01d0599556b8
# ╠═3653ee36-35a6-4e0a-8d46-4f8389381d45
# ╟─e774a20f-2d58-486a-ab71-6bde678b26f8
# ╠═589239fb-319c-40c2-af16-19025e7b28a2
# ╠═1f8f5bd0-20c8-4a52-9dac-4ceba18fcc06
# ╠═6fd88483-d005-4186-8dd2-82cea767ce90
# ╟─e3468c61-782b-4f55-a4a1-9d1883655d11
# ╠═bcda2043-f8c7-46bc-a5d4-b6f1f0883e9e
# ╠═f519626c-a3e8-4390-b3af-40b7beb665ed
# ╠═9a9b688c-94f0-4944-a9b2-21702073e0c7
# ╠═18d58341-0173-4eb1-9f01-cfa893088613
# ╟─941cd721-07d8-4a8f-9d75-42854e6e8edb
# ╠═4b763b58-862e-4c88-a7c9-fe0b1271c0b4
# ╠═2df82761-b9fe-4d37-b57c-1eabb0ffa8dd
# ╠═df46d106-f186-4900-9d3f-b711bc803707
# ╠═ab058d99-ce5f-4ed3-97bd-a62d2f258773
# ╠═13523326-a5f2-480d-9961-d23cd51192b8
# ╠═169197fe-983d-420b-8c56-353a65b28ddc
# ╟─4bad8b5c-e8b9-4ceb-97f4-41b4401d4f63
# ╠═ccabf5d2-5739-4284-a972-23c02a263a5c
# ╠═4e4cb513-1e88-4414-aa4d-a14d934874ce
# ╠═22b57aad-e886-4d36-bab8-baef5f3fabe6
# ╠═a4517d69-76e6-462a-9449-b31d80e34a8f
# ╠═ad5b07e5-75d0-4e03-a5d6-9ce4f1efd949
# ╟─0adc81ea-8678-42c2-a8b6-45fe4d26f4c4
# ╠═e6e1ea18-216a-41ae-8a1a-590793fcb669
# ╟─e98dee2e-a369-448e-bfe4-8fea0f318fa8
# ╠═3ca393d6-01c0-4f77-88ff-7c4f6388670e
# ╠═2768623f-7904-4674-a2ee-ad809cdd508b
# ╟─793c4d08-e2ee-4c9d-b7a0-11eaaddba895
# ╠═684c026a-b5d0-4694-8d29-a44b7cb0fd6c
# ╟─eeb3da97-72d5-4317-acb9-d28637a06d67
# ╠═a8d1c3e6-c020-495f-a443-07203b7dcd50
# ╠═b1b0690a-a1eb-11eb-1590-396d92c80c23
# ╟─d07e895f-c43b-46c4-be47-5af44bfda47a
# ╠═e8478e36-10fc-4e95-bf09-e217cad0cb15
# ╟─03af71ac-673b-459b-a931-a600b13d7ee6
