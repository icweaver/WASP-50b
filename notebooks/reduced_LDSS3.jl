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

# ╔═╡ e06eedec-5a47-4cd1-9f76-272c4c0476f2
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
end

# ╔═╡ 34ef4580-bb95-11eb-34c1-25893217f422
md"""
# Reduced

In this notebook we will examine the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from LDSS3.

$(Pl.TableOfContents(depth=4))
"""

# ╔═╡ a8e3b170-3fc2-4d12-b117-07bd37d27710
md"""
## Data extraction 🔳

The data from this instrument is stored in an `npy` file that contains the following fields:

* temporal
    * Table of external params (``N`` rows, where ``N`` is the number of timeseries points)
* squares
    * Dict with keys: `sky`, `shift`, `width`, `centroid`, and `stretch`. Each key contains a dict of N measurements, 1 for each aperture
* stellar
    * Dict of the aperture names, e.g., `aperture_830_689`
* spectral
    * Dict with two entries:
        * `wavelength`: Array of wavelengths in final grid, e.g., 5000-9799, with length ``W``
        * `fractionofapixel`: Amount of each pixel per wavelength bin
* comparisons
    * Similar to `stellar`, but a list of the comparison star square names only
* target
    * Similar to `stellar`, but a list of the target star square name(s) only
* cubes
    * Dict with keys: `raw_counts`, `peak`, `ok`, `sky`, `width`, and `centroid`. Each key contains an ``N \times W`` matrix of measurements, 1 matrix for each square.

!!! note
	We will use `wavelength` and `raw_counts` for our analysis here, where `raw_counts` is the sky subtracted, cross-dispersion integrated flux from each star.

Each cube (`LC`) can be selected from the following drop-down menu, and will be used for the rest of this analysis:
$(@bind dirpath Pl.Select(sort(glob("data/reduced/LDSS3/*/LC*.npy"))))
"""

# ╔═╡ 698ae5b0-7cd3-4055-a13e-e9aa3704ca12
md"""
We next define the mapping between the standard aperture names defined in `stellar` and the target/comparison star names used:
"""

# ╔═╡ d8236985-9a36-4357-ac78-7eb39dd0f080
names = OrderedDict(
	"aperture_324_803" => "WASP50",
	"aperture_153_1117" => "c06",
	"aperture_830_689" => "c15",
	"aperture_28_1189" => "c21",
)

# ╔═╡ b4f4e65b-bd50-4ee2-945f-7c130db21fdf
md"""
!!! note
	These stars were cross-matched with VizieR

	**Add to paper:** "This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: 10.26093/cds/vizier). The original description of the VizieR service was published in A&AS 143, 23"
"""

# ╔═╡ 0d017105-e3a3-4a0e-9e5f-be4200b8cc2c
CSV.File("data/reduced/LDSS3/asu.tsv", comment="#")

# ╔═╡ a687fa5d-1d98-4cb4-ae5d-978594d205dd
md"""
With this mapping, we then extract the `raw_counts` for each target and store them in `fluxes`:
"""

# ╔═╡ ae81cf2b-e2cf-42af-9bba-955155c63647
md"""
!!! note
	We convert from Float32 to Float64 to match the format of the IMACS data
"""

# ╔═╡ 839b83b2-97cd-4643-a4ce-9a03f3594b3a
md"""
#### Common flux

Next we extract the common wavelength grid (`wav`), along with the associated target `f_target` and comparison star flux (`f_comps`):
"""

# ╔═╡ 68610805-92fa-45a0-b5f6-7f04cb209c04
# Returns common column index range where rows of `A` != 0
function get_range(A)
	z = zero(eltype(A))
	rows = eachrow(A)
	first = maximum(findfirst.(!=(z), rows))
	last = minimum(findlast.(!=(z), rows))
	return first:last
end

# ╔═╡ 299dda0e-a214-45ca-9a68-947f60fcf404
md"""
## Stellar spectra 🌟

With the flux extracted for each object, we now turn to analyzing the resulting stellar spectra:
"""

# ╔═╡ a6088ea2-904f-4909-b1be-9470e7ec2010
med_std(A; dims=1) = (median(A, dims=dims), std(A, dims=dims)) .|> vec

# ╔═╡ 7d68ad39-3e39-48fa-939a-e56c6659d2b3
function spec_plot(ax, wav, A; norm=1.0, color=:blue, label="")
	μ, σ = med_std(A) ./ norm
	band(ax, wav, μ .- σ, μ .+ σ, color=(color, 0.25))
	lines!(ax, wav, μ, color=color, label=label)
end

# ╔═╡ bd937d51-17e9-4de3-a5d0-4c436d413940
md"""
## White light curves ⚪

Next, we will extract the integrated white light curves from these spectra. We integrate along the same wavelength binns used in the IMACS analysis:
"""

# ╔═╡ 9697e26b-b6d9-413b-869f-47bc2ab99919
wbins = readdlm("data/reduced/LDSS3/w50_bins_LDSS3.dat", comments=true)

# ╔═╡ cb805821-5d2e-484c-93a5-10a897d2cfe7
@bind window_width Pl.Slider(3:2:21, show_value=true)

# ╔═╡ d5c6d058-17c6-4cf0-97b8-d863b1529161
md"""
!!! note
	We divide the target WLC by each comparison star to i)minimize common systematics (e.g., air mass, local refractive atmospheric effects), and ii) make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ 89256633-3c12-4b94-b245-4fdda44d686c
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

# ╔═╡ f80347e8-dc5a-4b0c-a6c0-db5c12eadcbb
# Filter specified WLCs and return superset points
function filt_idxs(f_div_wlc, window_width; ferr=0.002)
	ntimes, ncomps = size(f_div_wlc)
	f_filt, diff = filt(f_div_wlc, window_width)
	bad_idxs = ∪(findall.(>(ferr), eachcol(diff))...) |> sort;
	use_idxs = deleteat!(collect(1:size(f_div_wlc, 1)), bad_idxs)
	return f_filt, use_idxs, bad_idxs
end

# ╔═╡ 9372c69a-0aad-4e6e-9ea3-e934fa09b758
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ 5fccf396-cd21-410f-9206-db93648b22f2
md"""
### Helper functions
"""

# ╔═╡ 8ec3c120-96e5-4fb4-b0c1-d8aae89b092f
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

# ╔═╡ 01e27cf6-0a1b-4741-828a-ef8bf7037ae9
LC = load_npz(dirpath, allow_pickle=true)

# ╔═╡ 9341c427-642a-4782-925a-6e07b91277a0
fluxes = Dict(
	names[k] => convert(Matrix{Float64}, v)
	for (k, v) ∈ LC["cubes"]["raw_counts"]
)

# ╔═╡ 651b5782-c1ae-45f0-8de1-079e13387ca9
# Non-zero wavelength ranges for each star
idx_lims = [get_range(f) for (_, f) in fluxes]

# ╔═╡ 22377e04-5a76-4d70-acf4-25d60f3b64a7
# Use the intersection for common wavelength range
common_wav_idxs = ∩(idx_lims...)

# ╔═╡ 9cdc7e8b-be59-44f9-9bee-4591e2ad788d
f_target = fluxes["WASP50"][:, common_wav_idxs]

# ╔═╡ 0d5749ac-4c94-4228-b721-83aaf84891c2
f_comps = cat(
	[
		fluxes[comp_name][:, common_wav_idxs]
		for comp_name in ["c06", "c15", "c21"]
	]...,
	dims=3
)

# ╔═╡ dc62887d-746b-4503-8547-7a6814de66a8
wav = LC["spectral"]["wavelength"][common_wav_idxs]

# ╔═╡ 8c82a78d-7382-40d4-a76b-02e7cd061d67
# Corresponding wav indxs
binned_wav_idxs = [
	findfirst(==(wbin[1]), wav):findfirst(==(wbin[2]), wav)
	for wbin in eachrow(wbins)
]

# ╔═╡ cf38810c-9b0c-4194-bea3-e0aa26e7cf98
binned_wav_idxs_range = binned_wav_idxs[begin][begin]:binned_wav_idxs[end][end]

# ╔═╡ 13385b21-fbd7-484d-a1ac-0687834f92c7
begin
	f_target_wlc = sum(f_target[:, binned_wav_idxs_range], dims=2)
	f_comps_wlc = sum(
		f_comps[:, binned_wav_idxs_range, :], dims=2
	) |> x -> dropdims(x, dims=2)
	f_div_wlc = f_target_wlc ./ f_comps_wlc
	f_div_WLC_norm = f_div_wlc ./ median(f_div_wlc, dims=1)
end

# ╔═╡ 470514e7-0f08-44a3-8519-5d704ea6b8d4
_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm, window_width);

# ╔═╡ 8b12f760-f294-4212-9b4e-88a886d84156
md"""
## Plot configs
"""

# ╔═╡ 21cbb4b9-1887-4575-8f4f-5e32b8404c6a
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

# ╔═╡ 45418bd3-74a3-4758-9fce-adddbeeec076
let
	fig = Figure(resolution=FIG_WIDE)

	fluxes = [f_target, [f_comps[:, :, i] for i in 1:3]...]
	labels = names.vals
	f_med = median(f_target)

	k = 1
	for j in 1:2, i in 1:2
		spec_plot(fig[i, j], wav, fluxes[k];
			norm = f_med,
			color = (COLORS[2]),
			label=labels[k],
		)
		k += 1
	end

	axs = reshape(copy(fig.content), 2, 2)
	axislegend.(axs)
	linkxaxes!(axs...)
	linkyaxes!(axs[begin, :]...)
	linkyaxes!(axs[end, :]...)
	hidexdecorations!.(axs[begin, :], grid=false)
	hideydecorations!.(axs[:, end], grid=false)

	fig[1:2, 0] = Label(fig, "Normalized counts", rotation=π/2)
	fig[end+1, 2:3] = Label(fig, "Index")
   #textsize = 30, font = noto_sans_bold, color = (:black, 0.25))

	fig |> Pl.as_svg
end

# ╔═╡ 798880fa-1e52-4758-a7f7-d3c6adec244a
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

# ╔═╡ 4b2ed9db-0a17-4e52-a04d-3a6a5bf2c054
let
	fig = Figure(resolution=FIG_TALL)
	
	comp_names = names.vals[begin+1:end]
	ncomps = length(comp_names)
	use_comps = comp_names # use all comps
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names))
	
	axs = [Axis(fig[i, 1]) for i in 1:ncomps]
	axs = reshape(copy(fig.content), ncomps, 1)
	
	plot_div_WLCS!(
		axs, f_div_WLC_norm, window_width, comp_names, use_comps_idxs
	)
	
	hidexdecorations!.(axs[begin:end-1], grid=false)
	linkaxes!(axs...)
	ylims!(axs[end], 0.97, 1.02)
	
	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2)
	fig[end+1, 1:end] = Label(fig, "Index")
	
	fig |> Pl.as_svg
end

# ╔═╡ f788835c-8e81-4afe-805e-4caf2d5e5d5b
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╟─34ef4580-bb95-11eb-34c1-25893217f422
# ╟─a8e3b170-3fc2-4d12-b117-07bd37d27710
# ╠═01e27cf6-0a1b-4741-828a-ef8bf7037ae9
# ╟─698ae5b0-7cd3-4055-a13e-e9aa3704ca12
# ╠═d8236985-9a36-4357-ac78-7eb39dd0f080
# ╟─b4f4e65b-bd50-4ee2-945f-7c130db21fdf
# ╠═0d017105-e3a3-4a0e-9e5f-be4200b8cc2c
# ╟─a687fa5d-1d98-4cb4-ae5d-978594d205dd
# ╠═9341c427-642a-4782-925a-6e07b91277a0
# ╟─ae81cf2b-e2cf-42af-9bba-955155c63647
# ╟─839b83b2-97cd-4643-a4ce-9a03f3594b3a
# ╠═651b5782-c1ae-45f0-8de1-079e13387ca9
# ╠═22377e04-5a76-4d70-acf4-25d60f3b64a7
# ╠═dc62887d-746b-4503-8547-7a6814de66a8
# ╠═68610805-92fa-45a0-b5f6-7f04cb209c04
# ╠═9cdc7e8b-be59-44f9-9bee-4591e2ad788d
# ╠═0d5749ac-4c94-4228-b721-83aaf84891c2
# ╟─299dda0e-a214-45ca-9a68-947f60fcf404
# ╠═45418bd3-74a3-4758-9fce-adddbeeec076
# ╠═7d68ad39-3e39-48fa-939a-e56c6659d2b3
# ╠═a6088ea2-904f-4909-b1be-9470e7ec2010
# ╟─bd937d51-17e9-4de3-a5d0-4c436d413940
# ╠═9697e26b-b6d9-413b-869f-47bc2ab99919
# ╠═8c82a78d-7382-40d4-a76b-02e7cd061d67
# ╠═cf38810c-9b0c-4194-bea3-e0aa26e7cf98
# ╠═13385b21-fbd7-484d-a1ac-0687834f92c7
# ╠═cb805821-5d2e-484c-93a5-10a897d2cfe7
# ╠═4b2ed9db-0a17-4e52-a04d-3a6a5bf2c054
# ╟─d5c6d058-17c6-4cf0-97b8-d863b1529161
# ╠═798880fa-1e52-4758-a7f7-d3c6adec244a
# ╠═f80347e8-dc5a-4b0c-a6c0-db5c12eadcbb
# ╠═89256633-3c12-4b94-b245-4fdda44d686c
# ╠═9372c69a-0aad-4e6e-9ea3-e934fa09b758
# ╠═470514e7-0f08-44a3-8519-5d704ea6b8d4
# ╟─5fccf396-cd21-410f-9206-db93648b22f2
# ╠═8ec3c120-96e5-4fb4-b0c1-d8aae89b092f
# ╟─8b12f760-f294-4212-9b4e-88a886d84156
# ╠═21cbb4b9-1887-4575-8f4f-5e32b8404c6a
# ╟─f788835c-8e81-4afe-805e-4caf2d5e5d5b
# ╠═e06eedec-5a47-4cd1-9f76-272c4c0476f2
