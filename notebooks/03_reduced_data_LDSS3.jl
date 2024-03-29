### A Pluto.jl notebook ###
# v0.19.8

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

# ╔═╡ f883b759-65fc-466e-9c8f-e4f941def935
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

# ╔═╡ 34ef4580-bb95-11eb-34c1-25893217f422
@mdx """
# Reduced data -- LDSS3C

In this notebook we will examine the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from LDSS3C.

!!! tip "Data download"
	
	```
	rclone sync -P ACCESS_box:WASP-50b/data/reduced data/reduced
	
	```
	
	* [Direct link](https://app.box.com/s/esq7gbpd7id98vzum1qub7twd4lba6zq)

	Outline:
	
	```
	reduced/
	├── [4.0K]  IMACS/
	│   └── [4.0K]  out_b/
	│       ├── [4.0K]  ut131219/
	│       │   └── [172M]  w50_131219_IMACS.pkl
	│       ├── [4.0K]  ut150927/
	│       │   └── [245M]  w50_150927_IMACS.pkl
	│       └── [4.0K]  ut161211/
	│           └── [199M]  w50_161211_IMACS.pkl
	├── [4.0K]  LDSS3C/
	│   ├── [3.2K]  catalog.tsv
	│   ├── [246K]  genericfinderchart.pdf
	│   ├── [4.0K]  LDSS3C_template/
	│   │   └── [4.0K]  WASP50/
	│   │       ├── [4.0K]  w50_150927_LDSS3_flat/
	│   │       │   ├── [ 27K]  eparams.dat
	│   │       │   ├── [4.0K]  wavelength/
	│   │       │   │   ├── [4.0K]  wbin0/
	│   │       │   │   │   ├── [ 17K]  comps.dat
	│   │       │   │   │   └── [ 13K]  lc.dat
	│   │       │   │   ├── [4.0K]  ⋅ ⋅ ⋅
	│   │       │   │   └── [4.0K]  wbin19/
	│   │       │   │       ├── [ 17K]  comps.dat
	│   │       │   │       └── [ 13K]  lc.dat
	│   │       │   └── [4.0K]  white-light/
	│   │       │       ├── [ 17K]  comps.dat
	│   │       │       └── [ 13K]  lc.dat
	│   │       └── [4.0K]  w50_150927_sp_LDSS3_flat/
	│   │           ├── [ 26K]  eparams.dat
	│   │           ├── [4.0K]  wavelength/
	│   │           │   ├── [4.0K]  wbin0/
	│   │           │   │   ├── [ 17K]  comps.dat
	│   │           │   │   └── [ 13K]  lc.dat
	│   │           │   ├── [4.0K]  wbin1/
	│   │           │   │   ├── [ 17K]  comps.dat
	│   │           │   │   └── [ 13K]  lc.dat
	│   │           │   ├── [4.0K]  ⋅ ⋅ ⋅
	│   │           │   └── [4.0K]  wbin14/
	│   │           │       ├── [ 17K]  comps.dat
	│   │           │       └── [ 13K]  lc.dat
	│   │           └── [4.0K]  white-light/
	│   │               ├── [ 17K]  comps.dat
	│   │               └── [ 13K]  lc.dat
	│   ├── [ 104]  targets.tsv
	│   ├── [4.0K]  ut150927_flat/
	│   │   ├── [122M]  LCs_ut150927.npy
	│   │   └── [412K]  specshifts.npy
	│   ├── [ 53K]  wasp50ls_mask.jpeg
	│   └── [1.1K]  wasp50ls.SMF
	└── [4.0K]  wbins/
	    ├── [ 288]  w50_bins.dat
	    ├── [ 302]  w50_bins_LDSS3.dat
	    ├── [ 225]  w50_bins_species.dat
	    └── [ 344]  w50_bins_ut131219.dat
	```
"""

# ╔═╡ e1294e56-8103-443b-80e9-118cac36750f
begin
	const DATA_DIR = "data/reduced/LDSS3C"
	const FIG_DIR = "figures/reduced"
	TableOfContents()
end

# ╔═╡ a8e3b170-3fc2-4d12-b117-07bd37d27710
@mdx """
## Data extraction 🔳

The data from this instrument are stored in an `npy` file that contains the following fields:

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
"""

# ╔═╡ 48713a10-6ba7-4de5-a147-85a9cb136dd8
@bind FPATH Select(glob("$(DATA_DIR)/ut150927_*/LCs_*.npy"))

# ╔═╡ 698ae5b0-7cd3-4055-a13e-e9aa3704ca12
@mdx """
We next define the mapping between the standard aperture names defined in `stellar` and the target/comparison star names used:
"""

# ╔═╡ d8236985-9a36-4357-ac78-7eb39dd0f080
obj_names = OrderedDict(
	"aperture_324_803" => "WASP50",
	"aperture_153_1117" => "c06",
	"aperture_830_689" => "c15",
	"aperture_28_1189" => "c21",
	"aperture_157_1116" => "c06",
	"aperture_830_693" => "c15",
	"aperture_35_1188" => "c21",
)

# ╔═╡ b4f4e65b-bd50-4ee2-945f-7c130db21fdf
@mdx """
!!! note
	These stars were cross-matched with VizieR to verify the aperture name assignments

	**Add to paper:** "This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: 10.26093/cds/vizier). The original description of the VizieR service was published in A&AS 143, 23"
"""

# ╔═╡ a687fa5d-1d98-4cb4-ae5d-978594d205dd
@mdx """
With this mapping, we then extract the `raw_counts` for each target and store them in `fluxes`:
"""

# ╔═╡ ae81cf2b-e2cf-42af-9bba-955155c63647
@mdx """
!!! note
	We convert from Float32 to Float64 to match the format of the IMACS data
"""

# ╔═╡ 839b83b2-97cd-4643-a4ce-9a03f3594b3a
@mdx """
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
@mdx """
## Stellar spectra ⭐

With the flux extracted for each object, we now turn to analyzing the resulting stellar spectra, selected from the wavelength bins scheme below:
"""

# ╔═╡ 9697e26b-b6d9-413b-869f-47bc2ab99919
@bind FPATH_WBINS let
	dirpath = "$(dirname(DATA_DIR))/wbins"
	Select(["$(dirpath)/w50_bins_LDSS3.dat", "$(dirpath)/w50_bins_species.dat"])
end

# ╔═╡ db933939-b8df-48b2-b53e-1dbd7ec1b07c
wbins = readdlm(FPATH_WBINS, comments=true)

# ╔═╡ 25903146-2f23-4e54-bd7d-aed1025f2fa5
fname_suff = let
	suff = "LDSS3C_" * basename(dirname(FPATH))
	occursin("species", FPATH_WBINS) ? (suff *= "_species") : suff
end

# ╔═╡ a6088ea2-904f-4909-b1be-9470e7ec2010
med_std(A; dims=1) = (median(A, dims=dims), std(A, dims=dims)) .|> vec

# ╔═╡ b9460e17-4034-4132-a57e-94721fa15d15
function spec_plot!(ax, wav, A; color=:blue, norm=1.0, label="")
	μ, σ = med_std(A) ./ norm
	band!(ax, wav, μ .- σ, μ .+ σ, color=(color, 0.25))
	lines!(ax, wav, μ; color, label)
	return μ, σ
end

# ╔═╡ bd937d51-17e9-4de3-a5d0-4c436d413940
@mdx """
## White-light curves 🌅

Next, we will extract the integrated white-light curves from these spectra. We integrate over the same wavelength bins used in the IMACS analysis:
"""

# ╔═╡ cb805821-5d2e-484c-93a5-10a897d2cfe7
@bind window_width PlutoUI.Slider(3:2:21, default=15, show_value=true)

# ╔═╡ 83263daf-a902-4414-850b-aa6949752fbb
cNames = obj_names.vals[2:4]

# ╔═╡ cd10cbf3-22f4-46cd-8345-cec3d141e3ca
use_comps = cNames

# ╔═╡ d386b581-a885-4461-9087-055269e77005
cNames_global = ("c06", "c13", "c15", "c18", "c20", "c21", "c23", "c28")

# ╔═╡ 274c809f-62f5-4470-a561-f005878e6fd2
function filt_curve(x; window_width=15, n_σ=2.0)
	x_med = mapwindow(median, x, window_width; border="reflect")
	x_err = mapwindow(std, x, window_width; border="reflect")
	x_diff = x - x_med
	x_diff_abs = abs.(x_diff)
	#bad_idxs = findall(x_diff .≥ (n_σ .* median(x_err)))
	bad_idxs = findall(x_diff_abs .≥ 0.002)
	return (; x_med, x_err, x_diff, x_diff_abs, bad_idxs)
end

# ╔═╡ 31732bf1-8797-497f-8669-40519898255e
mid_transit_times = Dict(
	"LDSS3C_ut150927_flat" => "2015-09-27 06:37",
	"LDSS3C_ut150927_flat_species" => "2015-09-27 06:37",
)

# ╔═╡ a82de6b2-35a6-4f90-8b38-1e50ffb1c586
date_fmt = dateformat"y-m-d H:M"

# ╔═╡ 6ce2ec33-bc6c-4f13-9904-2321658f9407
function compute_t_rel(t_py)
	t = pyconvert(Vector, t_py)
	t₀_utc = mid_transit_times[fname_suff]
	t₀ = DateTime(t₀_utc, date_fmt) |> datetime2julian
	return @. (t - t₀) * 24.0
end

# ╔═╡ 9372c69a-0aad-4e6e-9ea3-e934fa09b758
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ 4f71ba8d-bfa0-4adc-8d82-cd3bca8b6c14
use_comps_idxs = get_idx.(use_comps, Ref(cNames))

# ╔═╡ d5c6d058-17c6-4cf0-97b8-d863b1529161
@mdx """
!!! note
	We divide the target WLC by each comparison star to minimize common systematics (e.g., air mass, local refractive atmospheric effects), and to make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ 89256633-3c12-4b94-b245-4fdda44d686c
function filt(f_div_wlc, window_width; func=median, border="reflect")
	# Filtered light curves
	f_filt = mapslices(
		x -> mapwindow(func, x, window_width, border=border),
		f_div_wlc,
		dims = 1,
	)

	# Residualsa
	Δf = f_div_wlc - f_filt

	return f_filt, abs.(Δf), Δf
end

# ╔═╡ f80347e8-dc5a-4b0c-a6c0-db5c12eadcbb
# Filter specified WLCs and return superset points
function filt_idxs(f_div_wlc, window_width; ferr=0.002)
	ntimes, ncomps = size(f_div_wlc)
	f_filt, f_diff, _ = filt(f_div_wlc, window_width)
	bad_idxs = ∪(findall.(>(ferr), eachcol(f_diff))...) |> sort;
	use_idxs = deleteat!(collect(1:size(f_div_wlc, 1)), bad_idxs)
	return f_filt, use_idxs, bad_idxs
end

# ╔═╡ 97191c55-f673-46b4-82fd-147e7d150623
@mdx """
### Raw flux

Just as a quick check:
"""

# ╔═╡ 56efa000-a943-4256-94f4-f0ae4764f634
# let
# 	fig = Figure()
# 	ax = Axis(fig[1, 1];
# 		xlabel="Index",
# 		ylabel="Relative flux",
# 	)

# 	cNames = obj_names.vals[2:4]
# 	f_wlc_targ = f_target_wlc[:]
# 	f_wlc_comps = f_comps_wlc ./ mean(f_comps_wlc, dims=1)

# 	for (cName, col) in zip(sort(cNames), eachcol(f_wlc_comps))
# 		lines!(ax, col; label=cName)
# 	end

# 	lines!(ax, f_wlc_targ ./ mean(f_wlc_targ, dims=1);
# 		linewidth = 5,
# 		color = :darkgrey,
# 		label = "WASP-50",
# 	)

# 	axislegend()


# 	fig
# end

# ╔═╡ a5e742e5-fcca-40d7-b342-c6112e6899e5
@mdx """
## Binned light curves 🌈

We next integrate the target and comparison star flux within each bin defined above to build the binned light curves. We store these in `oLCw` and `cLCw`, where `oLCw` is an `ntimes` ``\times`` `nbins` matrix that holds the binned target flux, where `ntimes` is the number of timeseries points ``N``. Similarly `cLCw` is an `ntimes` ``\times`` `ncomps` ``\times`` `nbins` matrix that holds the comparison star flux, where `ncomps` is the number of comparison stars:
"""

# ╔═╡ 90b26a58-834a-445b-9242-d7b2e04be614
# For columnar data `A`, returns the sum of each row (for default `dims=2`)
# along the range of columns specified by `idx_range`
sum_flux(A, idx_range, dims=2) = sum(view(A, :, idx_range), dims=dims)[:, 1]

# ╔═╡ 5837dc0e-3537-4a90-bd2c-494e7b3e6bb7
@mdx """
With `oLCw` and `cLCw` now computed, we next compute `f_norm_w`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\times`` `ncomps` ``\times`` `nbins`, where for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve:
"""

# ╔═╡ af07dc54-eeb5-4fbe-8dd0-289dea97502a
@mdx """
## `GPTransmissionSpectra` inputs 🔩

Finally, we export the external parameters and light curves (in magnitude space) to be used for detrending in `GPTransmissionSpectra`:
"""

# ╔═╡ 88cc640d-b58d-4cde-b793-6c66e74f6b3a
@mdx """
#### External parameters
"""

# ╔═╡ 301ff07c-8dd5-403a-bae8-a4c38deeb331
target_name = "aperture_324_803" #"aperture_325_803" # "aperture_324_803"

# ╔═╡ dd8cf833-bb6f-4c70-a7c0-54fc7f08fdfc
rround(x) = round(x; digits=2)

# ╔═╡ d14ab9de-23b4-4647-a823-9b318bb734e9
median_eparam(param, cube, obj_name, wav_idxs) = cube[param][obj_name][:, wav_idxs] |>
	x -> median(x, dims=2) |> vec |> x -> convert(Vector{Float64}, x)

# ╔═╡ e6ed8da6-189a-4931-aac5-a4e8d0291723
fmt_float(x) = @sprintf "%.10f" x

# ╔═╡ 0be35b52-caea-4000-8cf8-ab99205bdb97
function template_dir(fpath)
	base_dir = "$(dirname(dirname(fpath)))/LDSS3C_template/WASP50"
	date, flat_status = split(basename(dirname(fpath)), '_')
	if occursin("species", fname_suff)
		f = "$(base_dir)/w50_$(date[3:end])_sp_LDSS3_$(flat_status)"
	else
		f = "$(base_dir)/w50_$(date[3:end])_LDSS3_$(flat_status)"
	end
	return f
end

# ╔═╡ 4af35a95-99b3-4186-a47d-169b9cbf927b
tdir = "$(template_dir(FPATH))"

# ╔═╡ 2341920d-6db8-4f08-b8ed-4907ceab7357
@mdx """
#### WLCs (magnitude space)
"""

# ╔═╡ b5ef9d95-1a2a-4828-8e27-89727f2e288b
@mdx """
#### Binned LCs (magnitude space)
"""

# ╔═╡ 123d0c63-f05a-4a7d-be16-6a3b9abac044
function f_to_med_mag(f)
	mag = -2.51 * log10.(f)
	return mag .- median(mag)
end

# ╔═╡ f788835c-8e81-4afe-805e-4caf2d5e5d5b
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ 36c1aa6d-cde9-4ff0-b55c-13f43e94256d
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig, pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 2d34e125-548e-41bb-a530-ba212c0ca17c
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
	load_npz(s; allow_pickle=false) = @pyeval("load_npz")(s, allow_pickle=allow_pickle)
	load_pickle(s) = @pyeval("load_pickle")(s)
end;

# ╔═╡ 01e27cf6-0a1b-4741-828a-ef8bf7037ae9
LC = load_npz(FPATH, allow_pickle=true);

# ╔═╡ 2d59bec3-9e2f-4f5f-b3f3-0963d546f3a0
LC_cubes = PyDict{String, PyDict{String, Matrix}}(LC["cubes"])

# ╔═╡ 9341c427-642a-4782-925a-6e07b91277a0
fluxes = Dict(obj_names[k] => v for (k, v) ∈ LC_cubes["raw_counts"])

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

# ╔═╡ 51563b2d-6d0b-4f20-a82e-36b58c3bca8f
LC_cubes["raw_counts"] |> typeof

# ╔═╡ 2f375e07-1e3e-473a-9488-ec482a42ec3c
LC_spectral = PyDict{String, Vector}(LC["spectral"])

# ╔═╡ dc62887d-746b-4503-8547-7a6814de66a8
wav = LC_spectral["wavelength"][common_wav_idxs]

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

# ╔═╡ c65c2298-e3a3-4666-be9d-73ee43d94847
FWHM, Trace_Center, Sky_Flux = median_eparam.(
	["width", "peak", "sky"],
	Ref(LC_cubes),
	Ref(target_name),
	Ref(common_wav_idxs[binned_wav_idxs_range])
)

# ╔═╡ 861cb600-5a97-496c-9a4d-8f848654f214
begin
	nbins = length(binned_wav_idxs)
	ntimes, ncomps = size(f_div_wlc)
	oLCw = Matrix{Float64}(undef, ntimes, nbins)
	cLCw = Array{Float64, 3}(undef, ntimes, nbins, ncomps)

	for (j, wav_idxs) in enumerate(binned_wav_idxs)
		oLCw[:, j] .= sum_flux(f_target, wav_idxs)
		for k in 1:ncomps
			cLCw[:, j, k] .= sum_flux(f_comps[:, :, k], wav_idxs)
		end
	end

	cLCw = permutedims(cLCw, [1, 3, 2]) # Match convention used in tepspec
end;

# ╔═╡ 823a4d10-698a-43fb-a4bb-1448670909eb
oLCw

# ╔═╡ 66b637dd-4a7f-4589-9460-67057f7945fd
cLCw

# ╔═╡ 49d04cf9-2bec-4350-973e-880e376ab428
begin
	offs = reshape(range(0, 0.3, length=nbins), 1, :) # Arbitrary offsets for clarity
	f_norm_w = Array{Float64}(undef, ntimes, ncomps, nbins)
	for c_i in 1:ncomps
		f_w = oLCw ./ cLCw[:, c_i, :]
		f_norm_w[:, c_i, :] .= f_w ./ median(f_w, dims=1) .+ offs
	end
	baselines = ones(ntimes, nbins) .+ offs # Reference baselines
end;

# ╔═╡ 8bb36e39-ade4-4799-ab91-92caf38021b4
PyDict{String, Vector}(LC["spectral"])

# ╔═╡ c03cb527-d16d-47aa-ab63-6970f4ff0b1f
Times, Airmass = let
	vals = DataFrame(PyPandasDataFrame(LC["temporal"].to_pandas()))
	vals.bjd, vals.airmass
end

# ╔═╡ 8225c739-f825-43d7-9fd5-064506c619b5
t_rel = compute_t_rel(Times)

# ╔═╡ 52a7c738-6f53-46b7-a164-06d439f33d1c
begin
	n_σ = 2.0
	# Only apply filter to specified comp star divided WLCs
	filt_curves = filt_curve.(
		eachcol(f_div_WLC_norm[:, use_comps_idxs]);
		window_width,
		n_σ,
	)
	med_models = [x.x_med for x ∈ filt_curves]
	bad_idxs = [x.bad_idxs for x ∈ filt_curves]
	bad_idxs_common = ∪((x.bad_idxs for x ∈ filt_curves)...) |> sort
	use_idxs_common = deleteat!(collect(1:length(t_rel)), bad_idxs_common)
end;

# ╔═╡ ded63b4b-61b6-41b6-98d4-d13166bce76a
@with_terminal begin
	println(bad_idxs_common .- 1)
	println(length(bad_idxs_common))
end

# ╔═╡ bb186edb-bf2e-4984-9a35-d9f373d9d3d5
use_idxs = deleteat!(collect(1:size(f_div_wlc, 1)), bad_idxs_common)

# ╔═╡ bc54942e-38ef-4919-a3d4-28d5f4db8487
comps = let
	mag = -2.51 * log10.(f_comps_wlc)
	mag .-= median(mag, dims=1)
	mag[use_idxs, :]
end

# ╔═╡ 5110de9a-3721-4043-b8b7-493daacb4137
target_binned_mags = mapslices(f_to_med_mag, oLCw, dims=1)[use_idxs, :]

# ╔═╡ 53f5a645-93e0-499a-bb36-e4ff0153a63c
comp_binned_mags = mapslices(f_to_med_mag, cLCw, dims=1)[use_idxs, :, :]

# ╔═╡ f56c04c9-a684-4264-80ce-ba41d3f6e527
let 
	Z = Airmass
	b, m, e = Z[begin], Z[(begin+end) ÷ 2], Z[end]
	println("$(rround(b)) \\to $(rround(m)) \\to $(rround(e))")
end

# ╔═╡ 354580e4-0aa9-496f-b024-665025a2eeda
lc = let
	med_mag = f_to_med_mag(f_target_wlc |> vec)
	hcat(Times, med_mag, zeros(length(Times)))[use_idxs, :]
end

# ╔═╡ 898a390b-49f7-45f4-b1a1-b22922d69a29
let
    savepath = "$(tdir)/white-light"
	rm(savepath, force=true, recursive=true)
	mkpath(savepath)
	f = "$(savepath)/lc.dat"
	writedlm(f, lc, ",    ")
	@info "Saved to $(f)"
	f = "$(savepath)/comps.dat"
	writedlm(f, comps, ",    ")
	@info "Saved to $(f)"
end

# ╔═╡ 631c4d02-58cc-4c70-947f-22c8e8b11015
let
    savepath = "$(tdir)/wavelength"
	rm(savepath, force=true, recursive=true)
	for i in 1:nbins
		save_path_w = "$(savepath)/wbin$(i-1)"
		mkpath("$(save_path_w)")
		lc_w = hcat(
			Times[use_idxs], target_binned_mags[:, i], zeros(length(Times[use_idxs]))
		)
		f = "$(save_path_w)/lc.dat"
		writedlm(f, lc_w, ",    ")
		@info "Saved to $(f)"
		f = "$(save_path_w)/comps.dat"
		writedlm(f, comp_binned_mags[:, :, i], ",    ")
		@info "Saved to $(f)"
	end
end

# ╔═╡ 079c3915-33af-40db-a544-28453732c372
specshifts = load_npz(
	"$(dirname(FPATH))/specshifts.npy", allow_pickle=true
);

# ╔═╡ e4960d1a-8e33-478a-8100-d1838782938d
Delta_Wav = pyconvert(Dict{String, Float64}, specshifts["shift"][target_name]) |>
		sort |> values |> collect |> x -> convert(Vector{Float64}, x)

# ╔═╡ e4388fba-64ef-4588-a1ed-283da2f52196
df_eparams = DataFrame(
	(;Times, Airmass, Delta_Wav, FWHM, Sky_Flux, Trace_Center)
) |> df -> mapcols(col -> fmt_float.(col), df)[use_idxs, :]

# ╔═╡ 2aba612a-7599-4a2d-9ff0-2fd398c2a0db
let
    savepath = template_dir(FPATH)
	rm(savepath, force=true, recursive=true)
	mkpath(savepath)
	f = "$(tdir)/eparams.dat"
	CSV.write(f, df_eparams, delim=",    ")
	@info "Saved to $(f)"
end

# ╔═╡ c911cecd-0747-4cd1-826f-941f2f58091c
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = 72 .* (6, 8)
	const FIG_WIDE = 72 .* (12, 6)
	const FIG_LARGE = 72 .* (12, 12)
	const COLORS_SERIES = categorical_colors(:seaborn_colorblind, 9)
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
			Scatter = (linewidth=10, strokewidth=0),
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

# ╔═╡ 1ec8c68f-4e76-4764-a520-320aebbefc27
cName_color = OrderedDict(
	cName => color
	for (cName, color) ∈ zip(["WASP50"; cNames_global...], COLORS_SERIES)
);

# ╔═╡ 45418bd3-74a3-4758-9fce-adddbeeec076
let
	wbins_odd = @view wbins[begin:2:end, :]
	wbins_even = @view wbins[begin+1:2:end, :]
	# Larger font for two-column
	fig = Figure(resolution=FIG_WIDE, fontsize=24)
	ax = Axis(fig[1, 1];
		xlabel = "Wavelength Å",
		ylabel = "Relative flux",
		xlabelsize = 24,
		ylabelsize = 24,
		limits = (4000, 11_000, 0.0, 2.8),
		xticks = (4_500:1000:10_500),
		yticks = 0:0.5:3.0,
	)

	fluxes = [f_target, [f_comps[:, :, i] for i in 1:3]...]
	labels = obj_names.vals
	norm = 40_000.0 #median(f_target)

	vspan!(ax, wbins_odd[:, 1], wbins_odd[:, 2], color=(:black, 0.25))
	vspan!(ax, wbins_even[:, 1], wbins_even[:, 2], color=(:darkgrey, 0.25))
	
	for (i, (label, f)) in enumerate(zip(labels, fluxes))
		μ, σ = spec_plot!(ax, wav, f;
			color = cName_color[label],
			norm,
			label,
		)
	end
	axislegend("Transit 2 (LDSS3C)", halign=:right, gridshalign=:right)

	savefig(fig, "$(FIG_DIR)/extracted_spectra_$(fname_suff).pdf")

	fig
end

# ╔═╡ cbcb7cb7-4ada-4c5b-94ae-b33b779d580f
function plot_div_WLCS!(axs, t_rel, f; window_width, cNames, n_σ)
	k = 1
	color = :darkgrey
	z = 1
	for (i, cName) ∈ enumerate(cNames_global)
		if (z ≤ ncomps) && cNames_global[i] == cNames[z]
			# All points
			if cName ∈ ("c06", "c15", "c21") # LDSS3C comps
				c_text = COLORS[end]
			else
				c_text = :darkgrey
			end
			scatter!(axs[i], t_rel, f[:, z]; markersize=5, color=(color, 0.25))
			text!(axs[i], "$(cName)";
				position =(3, 0.98),
				align = (:right, :center),
				color = c_text,
				textsize = 24,
			)
			# Used points
			if cName ∈ use_comps
				# scatter!(axs[i], t_rel[bad_idxs[k]], f[bad_idxs[k], z];
				# 	marker = '⭘',
				# 	markersize = 20,
				# 	linewidth = 3.0,
				# 	color = COLORS_SERIES[4],
				# )
				scatter!(axs[i], t_rel[use_idxs_common], f[use_idxs_common, z];
					markersize = 5,
					color,
				)
				lines!(axs[i], t_rel, med_models[k];
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

	return bad_idxs_common
end

# ╔═╡ 1bb71101-adc8-4e9b-9354-d7411b131920
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

	plot_div_WLCS!(axs, t_rel, f_div_WLC_norm; window_width, cNames, n_σ)

	linkaxes!(axs...)
	hidexdecorations!.(axs[begin:end-1, :], grid=false)
	hideydecorations!.(axs[:, begin+1:end], grid=false)

	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2, textsize=24)
	fig[end+1, 2:end] = Label(fig, "Time from estimated mid-transit (hours)", textsize=24)

	Label(fig[0, end], "Transit 2 (LDSS3C)";
		tellwidth = false,
		halign = :right,
		textsize=24,
	)

	savefig(fig, "$(FIG_DIR)/div_wlcs_$(fname_suff).pdf")

	fig
end

# ╔═╡ 2419e060-f5ab-441b-9ec2-51ce4e57e319
function plot_BLCs(t, datas, models, wbins, errs, comp_name; offset=0.3)
	fig = Figure(resolution=FIG_LARGE)
	median_prec = round(Int, median(errs))

	ax_left = Axis(fig[1, 1], title = "Divided BLCs")
	ax_right = Axis(fig[1, 2], title = "Residuals")
	ax_label = Axis(fig[1, 3], title = "Median precision: $(median_prec) ppm")
	axs = reshape(copy(fig.content), 1, 3)
	linkaxes!(axs...)

	# Color palette
	N_bins = size(datas, 2)
	resids = datas - models
	colors = resample_cmap(:Spectral_4, N_bins) |> reverse

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
		scatter!(ax_left, t, data, strokewidth=0, markersize=5, color=color)
		lines!(ax_left, t, model, linewidth=3, color=0.75*color)

		scatter!(ax_right, t, baseline + resid, markersize=5, color=color)
		lines!(ax_right, t, baseline, linewidth=3, color=0.75*color)

		scatter!(ax_right, t, baseline + resid, markersize=5, color=color)
		lines!(ax_right, t, baseline, linewidth=3, color=0.75*color)
		text!(ax_label, "$(wbin[1]) - $(wbin[2]) Å, $(err) ppm";
			position = (-1.5, baseline[1]),
			textsize = 16,
			align = (:left, :center),
			offset = (0, 2),
			color = 0.75*color,
		)
	end

	hideydecorations!.(axs[:, 2:3], grid=false)
	hidespines!(axs[end])
	hidedecorations!(axs[end])
	ylims!(ax_left, 0.95, 1.34)

	fig[1:2, 0] = Label(fig, "Relative flux + offset", rotation=π/2)
	fig[end, 2:3] = Label(fig, "Time from estimated mid-transit (hours)")

	savefig(fig, "$(FIG_DIR)/div_blcs_$(fname_suff)_$(comp_name).pdf")

	fig
end

# ╔═╡ 65c91d63-10e7-41ab-9c21-f136f4c5cb96
begin
	blc_plots = OrderedDict()
	for comp_idx ∈ use_comps_idxs
		t = t_rel[use_idxs_common]
		datas = f_norm_w[:, comp_idx, :]
		cName = cNames[comp_idx]
		f_med, _, f_diff = filt(datas, window_width)
		p = plot_BLCs(
			t,
			datas[use_idxs, :],
			f_med[use_idxs, :],
			wbins,
			round.(Int, reshape(std(f_diff[use_idxs, : ], dims=1), :, 1) * 1e6),
			cName,
		)
		blc_plots[cName] = p
	end
end

# ╔═╡ 6dece5df-0e78-42f1-a557-0a5445350b28
@bind cName Select(blc_plots.keys)

# ╔═╡ ee489959-654a-4666-854d-79896832dbc7
@mdx """
target / $(cName)
"""

# ╔═╡ 7fa35566-d327-4319-9ae9-17c4c9825e05
blc_plots[cName]

# ╔═╡ Cell order:
# ╟─34ef4580-bb95-11eb-34c1-25893217f422
# ╠═e1294e56-8103-443b-80e9-118cac36750f
# ╟─a8e3b170-3fc2-4d12-b117-07bd37d27710
# ╟─48713a10-6ba7-4de5-a147-85a9cb136dd8
# ╠═01e27cf6-0a1b-4741-828a-ef8bf7037ae9
# ╟─698ae5b0-7cd3-4055-a13e-e9aa3704ca12
# ╠═d8236985-9a36-4357-ac78-7eb39dd0f080
# ╟─b4f4e65b-bd50-4ee2-945f-7c130db21fdf
# ╟─a687fa5d-1d98-4cb4-ae5d-978594d205dd
# ╠═9341c427-642a-4782-925a-6e07b91277a0
# ╠═51563b2d-6d0b-4f20-a82e-36b58c3bca8f
# ╠═2d59bec3-9e2f-4f5f-b3f3-0963d546f3a0
# ╟─ae81cf2b-e2cf-42af-9bba-955155c63647
# ╟─839b83b2-97cd-4643-a4ce-9a03f3594b3a
# ╠═651b5782-c1ae-45f0-8de1-079e13387ca9
# ╠═22377e04-5a76-4d70-acf4-25d60f3b64a7
# ╠═2f375e07-1e3e-473a-9488-ec482a42ec3c
# ╠═dc62887d-746b-4503-8547-7a6814de66a8
# ╠═8bb36e39-ade4-4799-ab91-92caf38021b4
# ╠═68610805-92fa-45a0-b5f6-7f04cb209c04
# ╠═9cdc7e8b-be59-44f9-9bee-4591e2ad788d
# ╠═0d5749ac-4c94-4228-b721-83aaf84891c2
# ╟─299dda0e-a214-45ca-9a68-947f60fcf404
# ╟─9697e26b-b6d9-413b-869f-47bc2ab99919
# ╟─db933939-b8df-48b2-b53e-1dbd7ec1b07c
# ╠═25903146-2f23-4e54-bd7d-aed1025f2fa5
# ╠═1ec8c68f-4e76-4764-a520-320aebbefc27
# ╠═45418bd3-74a3-4758-9fce-adddbeeec076
# ╠═b9460e17-4034-4132-a57e-94721fa15d15
# ╠═a6088ea2-904f-4909-b1be-9470e7ec2010
# ╟─bd937d51-17e9-4de3-a5d0-4c436d413940
# ╠═8c82a78d-7382-40d4-a76b-02e7cd061d67
# ╠═cf38810c-9b0c-4194-bea3-e0aa26e7cf98
# ╠═13385b21-fbd7-484d-a1ac-0687834f92c7
# ╠═cb805821-5d2e-484c-93a5-10a897d2cfe7
# ╠═ded63b4b-61b6-41b6-98d4-d13166bce76a
# ╠═cd10cbf3-22f4-46cd-8345-cec3d141e3ca
# ╠═83263daf-a902-4414-850b-aa6949752fbb
# ╠═4f71ba8d-bfa0-4adc-8d82-cd3bca8b6c14
# ╠═d386b581-a885-4461-9087-055269e77005
# ╠═1bb71101-adc8-4e9b-9354-d7411b131920
# ╠═8225c739-f825-43d7-9fd5-064506c619b5
# ╠═274c809f-62f5-4470-a561-f005878e6fd2
# ╠═31732bf1-8797-497f-8669-40519898255e
# ╠═a82de6b2-35a6-4f90-8b38-1e50ffb1c586
# ╠═6ce2ec33-bc6c-4f13-9904-2321658f9407
# ╠═9372c69a-0aad-4e6e-9ea3-e934fa09b758
# ╟─d5c6d058-17c6-4cf0-97b8-d863b1529161
# ╠═52a7c738-6f53-46b7-a164-06d439f33d1c
# ╠═cbcb7cb7-4ada-4c5b-94ae-b33b779d580f
# ╠═f80347e8-dc5a-4b0c-a6c0-db5c12eadcbb
# ╠═89256633-3c12-4b94-b245-4fdda44d686c
# ╟─97191c55-f673-46b4-82fd-147e7d150623
# ╠═56efa000-a943-4256-94f4-f0ae4764f634
# ╟─a5e742e5-fcca-40d7-b342-c6112e6899e5
# ╠═861cb600-5a97-496c-9a4d-8f848654f214
# ╠═90b26a58-834a-445b-9242-d7b2e04be614
# ╠═823a4d10-698a-43fb-a4bb-1448670909eb
# ╠═66b637dd-4a7f-4589-9460-67057f7945fd
# ╟─5837dc0e-3537-4a90-bd2c-494e7b3e6bb7
# ╠═49d04cf9-2bec-4350-973e-880e376ab428
# ╠═65c91d63-10e7-41ab-9c21-f136f4c5cb96
# ╟─6dece5df-0e78-42f1-a557-0a5445350b28
# ╟─ee489959-654a-4666-854d-79896832dbc7
# ╟─7fa35566-d327-4319-9ae9-17c4c9825e05
# ╠═2419e060-f5ab-441b-9ec2-51ce4e57e319
# ╟─af07dc54-eeb5-4fbe-8dd0-289dea97502a
# ╠═bb186edb-bf2e-4984-9a35-d9f373d9d3d5
# ╟─88cc640d-b58d-4cde-b793-6c66e74f6b3a
# ╠═301ff07c-8dd5-403a-bae8-a4c38deeb331
# ╠═c03cb527-d16d-47aa-ab63-6970f4ff0b1f
# ╠═dd8cf833-bb6f-4c70-a7c0-54fc7f08fdfc
# ╠═f56c04c9-a684-4264-80ce-ba41d3f6e527
# ╠═c65c2298-e3a3-4666-be9d-73ee43d94847
# ╠═d14ab9de-23b4-4647-a823-9b318bb734e9
# ╠═079c3915-33af-40db-a544-28453732c372
# ╠═e4960d1a-8e33-478a-8100-d1838782938d
# ╠═e4388fba-64ef-4588-a1ed-283da2f52196
# ╠═e6ed8da6-189a-4931-aac5-a4e8d0291723
# ╠═2aba612a-7599-4a2d-9ff0-2fd398c2a0db
# ╠═4af35a95-99b3-4186-a47d-169b9cbf927b
# ╠═0be35b52-caea-4000-8cf8-ab99205bdb97
# ╟─2341920d-6db8-4f08-b8ed-4907ceab7357
# ╠═354580e4-0aa9-496f-b024-665025a2eeda
# ╠═bc54942e-38ef-4919-a3d4-28d5f4db8487
# ╠═898a390b-49f7-45f4-b1a1-b22922d69a29
# ╟─b5ef9d95-1a2a-4828-8e27-89727f2e288b
# ╠═5110de9a-3721-4043-b8b7-493daacb4137
# ╠═53f5a645-93e0-499a-bb36-e4ff0153a63c
# ╠═631c4d02-58cc-4c70-947f-22c8e8b11015
# ╠═123d0c63-f05a-4a7d-be16-6a3b9abac044
# ╟─f788835c-8e81-4afe-805e-4caf2d5e5d5b
# ╟─36c1aa6d-cde9-4ff0-b55c-13f43e94256d
# ╠═2d34e125-548e-41bb-a530-ba212c0ca17c
# ╠═c911cecd-0747-4cd1-826f-941f2f58091c
# ╠═f883b759-65fc-466e-9c8f-e4f941def935
