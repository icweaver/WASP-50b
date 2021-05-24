### A Pluto.jl notebook ###
# v0.14.5

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

# ╔═╡ b1b0690a-a1eb-11eb-1590-396d92c80c23
begin
	import PlutoUI as Pl
	using CairoMakie
	using Glob
	using PyCall
	using Statistics
	using Colors
	using ImageFiltering
	using CSV
	using DataFrames
	using DelimitedFiles
	using OrderedCollections
	using Printf
end

# ╔═╡ ee24f7df-c4db-4065-afe9-10be80cbcd6b
md"""
# Reduced

In this notebook we will examine the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from IMACS.

$(Pl.TableOfContents(depth=4))
"""

# ╔═╡ 9d180c21-e634-4a1e-8430-bdd089262f66
md"""
## Data extraction 🔳

The main data product from our custom pipeline for this instrument is a pickle file with the following naming scheme: `LCs_<target>_<wavelength bin scheme>.pkl`. 

Each cube (`LC`) can be selected from the following drop-down menu, and will be used for the rest of this analysis:
$(@bind fpath Pl.Select(sort(glob("data/reduced/IMACS/*/*.pkl"))))
"""

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

# ╔═╡ bd2cdf33-0c41-4948-82ab-9a28929f72b3
LC = load_pickle(fpath)

# ╔═╡ e774a20f-2d58-486a-ab71-6bde678b26f8
md"""
## Stellar spectra 🌟

With the flux extracted for each object, we now turn to analyzing the resulting stellar spectra:
"""

# ╔═╡ 52af8a0f-7ebf-4542-86fb-fe6a1b3565cc
median(LC["spectra"]["WASP50"])

# ╔═╡ 6fd88483-d005-4186-8dd2-82cea767ce90
med_std(A; dims=1) = (median(A, dims=dims), std(A, dims=dims)) .|> vec

# ╔═╡ 1f8f5bd0-20c8-4a52-9dac-4ceba18fcc06
function spec_plot!(ax, wav, A; norm=1.0, label="")
	μ, σ = med_std(A) ./ norm
	band!(ax, wav, μ .- σ, μ .+ σ)
	lines!(ax, wav, μ, label=label)
end

# ╔═╡ 589239fb-319c-40c2-af16-19025e7b28a2
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	wav = LC["spectra"]["wavelengths"]
	f_norm = median(LC["spectra"]["WASP50"])
	
	for (name, f) in sort(LC["spectra"])
		if name != "wavelengths"
			spec_plot!(ax, wav, f, norm=f_norm, label=name)
		end
	end
	
	axislegend()
	
	xlims!(ax, 3_500, 11_000)
	ylims!(ax, 0, 2.2)
	
	fig |> Pl.as_svg
end

# ╔═╡ e3468c61-782b-4f55-a4a1-9d1883655d11
md"""
## White light curves ⚪

Next, we will extract the integrated white light curves from these spectra, integrated over the specified wavelength bins:
"""

# ╔═╡ 18d58341-0173-4eb1-9f01-cfa893088613
begin
	comp_names = LC["cNames"]
	sorted_cName_idxs = sortperm(comp_names)
	
	f_div_WLC = LC["oLC"] ./ LC["cLC"][:, sorted_cName_idxs]
	f_div_WLC_norm = f_div_WLC ./ median(f_div_WLC, dims=1)
end

# ╔═╡ 941cd721-07d8-4a8f-9d75-42854e6e8edb
md"""
!!! warning
	In general, the comparison stars names (`cNames`) are not stored in alphanumeric order by default. For convenience, we ensure this sorting with `sortperm`, so that the first column corresponds to the first comparison star, the second to the second comparison star, an so on.

We next plot these light curves and identified outliers below:
"""

# ╔═╡ ab058d99-ce5f-4ed3-97bd-a62d2f258773
@bind window_width Pl.Slider(3:2:21, show_value=true)

# ╔═╡ 4b763b58-862e-4c88-a7c9-fe0b1271c0b4
use_comps = ["c15", "c21"]

# ╔═╡ 169197fe-983d-420b-8c56-353a65b28ddc
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ 4bad8b5c-e8b9-4ceb-97f4-41b4401d4f63
md"""
!!! note
	We divide the target WLC by each comparison star to minimize common systematics (e.g., air mass, local refractive atmospheric effects), and to make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ dc044a72-4706-49e2-94a8-c828a6bf7de0
function filt(f_div_wlc, window_width; func=median, border="reflect")
	# Filtered light curves
	f_filt = mapslices(
		x -> mapwindow(func, x, window_width, border=border),
		f_div_wlc,
		dims=1,
	)
	
	# Residuals
	diff = abs.(f_div_wlc - f_filt)
	
	return f_filt, diff
end

# ╔═╡ a4517d69-76e6-462a-9449-b31d80e34a8f
# Filter specified WLCs and return superset points
function filt_idxs(f_div_wlc, window_width; ferr=0.002)
	ntimes, ncomps = size(f_div_wlc)
	f_filt, diff = filt(f_div_wlc, window_width)
	bad_idxs = ∪(findall.(>(ferr), eachcol(diff))...) |> sort;
	use_idxs = deleteat!(collect(1:size(f_div_wlc, 1)), bad_idxs)
	return f_filt, use_idxs, bad_idxs
end

# ╔═╡ df46d106-f186-4900-9d3f-b711bc803707
Pl.with_terminal() do
	use_comps = use_comps
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names))
	_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm[:, use_comps_idxs], window_width)
	# Because python
	println(bad_idxs .- 1)
	println(use_comps_idxs .- 1)
end

# ╔═╡ e98dee2e-a369-448e-bfe4-8fea0f318fa8
md"""
## Binned light curves 🌈
"""

# ╔═╡ a604260b-902e-44a1-88ea-440438e582ed
md"""
## Plot configs
"""

# ╔═╡ 5db4a2f2-1c0d-495a-8688-40fc9e0ccd02
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

# ╔═╡ ccabf5d2-5739-4284-a972-23c02a263a5c
function plot_div_WLCS!(
	axs, f_div_wlc, window_width, cNames, use_comp_idxs; ferr=0.002
)
	use_comps = cNames[use_comp_idxs]
	
	# Only apply filter to specified comp star divided WLCs
	f_filt, use_idxs, bad_idxs = filt_idxs(
		f_div_wlc[:, use_comp_idxs], window_width; ferr=ferr
	)
	
	idxs = 1:size(f_div_wlc, 1)
	c = COLORS[end]
	k = 1
	for (i, cName) ∈ enumerate(cNames)		
		# All points
		scatter!(axs[i], idxs, f_div_wlc[:, i];
			color = (c, 0.3),
			strokewidth = 0,
			label = "$cName",
		)
		# Used points
		if cName ∈ use_comps
			scatter!(axs[i], idxs[use_idxs], f_div_wlc[use_idxs, i];
				color = c,
				strokewidth = 0,
			)
			lines!(axs[i], idxs, f_filt[:, k];
				color = COLORS[end-2],
				linewidth = 2,
			)
			k += 1
		end
		
		axislegend(axs[i])
	end
end

# ╔═╡ 13523326-a5f2-480d-9961-d23cd51192b8
let
	fig = Figure(resolution=FIG_TALL)
	
	ncomps = length(comp_names)
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names))
	
	axs = [Axis(fig[i, j]) for i ∈ 1:4, j ∈ 1:2]
	axs = reshape(copy(fig.content), 4, 2)
	
	plot_div_WLCS!(
		axs, f_div_WLC_norm, window_width, comp_names, use_comps_idxs
	)
	
	linkaxes!(axs...)
	hidexdecorations!.(axs[begin:end-1, :], grid=false)
	hideydecorations!.(axs[:, begin+1:end], grid=false)
	ylims!(axs[end], 0.97, 1.02)
	
	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2)
	fig[end+1, 1:end] = Label(fig, "Index")
	
	fig |> Pl.as_svg
end

# ╔═╡ eeb3da97-72d5-4317-acb9-d28637a06d67
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╠═ee24f7df-c4db-4065-afe9-10be80cbcd6b
# ╟─9d180c21-e634-4a1e-8430-bdd089262f66
# ╠═bd2cdf33-0c41-4948-82ab-9a28929f72b3
# ╟─66052b03-35a0-4877-abef-f525766fd332
# ╟─7bfc971c-8737-49ad-adec-ac57d176f10e
# ╟─1c3e8cb3-2eff-47c2-8c17-01d0599556b8
# ╠═3653ee36-35a6-4e0a-8d46-4f8389381d45
# ╟─e774a20f-2d58-486a-ab71-6bde678b26f8
# ╠═589239fb-319c-40c2-af16-19025e7b28a2
# ╠═52af8a0f-7ebf-4542-86fb-fe6a1b3565cc
# ╠═1f8f5bd0-20c8-4a52-9dac-4ceba18fcc06
# ╠═6fd88483-d005-4186-8dd2-82cea767ce90
# ╟─e3468c61-782b-4f55-a4a1-9d1883655d11
# ╠═18d58341-0173-4eb1-9f01-cfa893088613
# ╟─941cd721-07d8-4a8f-9d75-42854e6e8edb
# ╠═ab058d99-ce5f-4ed3-97bd-a62d2f258773
# ╠═4b763b58-862e-4c88-a7c9-fe0b1271c0b4
# ╠═13523326-a5f2-480d-9961-d23cd51192b8
# ╠═169197fe-983d-420b-8c56-353a65b28ddc
# ╟─4bad8b5c-e8b9-4ceb-97f4-41b4401d4f63
# ╠═ccabf5d2-5739-4284-a972-23c02a263a5c
# ╠═a4517d69-76e6-462a-9449-b31d80e34a8f
# ╠═dc044a72-4706-49e2-94a8-c828a6bf7de0
# ╠═df46d106-f186-4900-9d3f-b711bc803707
# ╟─e98dee2e-a369-448e-bfe4-8fea0f318fa8
# ╟─a604260b-902e-44a1-88ea-440438e582ed
# ╠═5db4a2f2-1c0d-495a-8688-40fc9e0ccd02
# ╟─eeb3da97-72d5-4317-acb9-d28637a06d67
# ╠═b1b0690a-a1eb-11eb-1590-396d92c80c23
