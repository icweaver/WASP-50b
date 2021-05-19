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

# ╔═╡ 442946ba-9561-4223-8223-ba513ac7bc45
using Printf

# ╔═╡ b1b0690a-a1eb-11eb-1590-396d92c80c23
begin
	import PlutoUI as pl
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

# ╔═╡ ee24f7df-c4db-4065-afe9-10be80cbcd6b
md"""
# Reduced

In this notebook we will examine the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from IMACS and LDSS3.

$(pl.TableOfContents(depth=4))
"""

# ╔═╡ 9d180c21-e634-4a1e-8430-bdd089262f66
md"""
## Data extraction 🔳

We demonstrate the extraction process for the data products from each instrument below.
"""

# ╔═╡ 454cab9c-7ea0-4b0a-9622-d29d8ba0395b
md"""
### IMACS

The main data product from our custom pipeline for this instrument is a pickle file with the following naming scheme: `LCs_<target>_<wavelength bin scheme>.pkl`. 

Each cube can be selected from the following drop-down menu:
$(@bind data_path_IMACS pl.Select(sort(glob("data/reduced/IMACS/*/*.pkl"))))
"""

# ╔═╡ 66052b03-35a0-4877-abef-f525766fd332
md"""
!!! tip
	A description of each field can be found on our repo page here (PUBLIC LINK?)
"""

# ╔═╡ 470c738f-c2c8-4f56-b936-e08c43c161b1
md"""
### LDSS3

This data is stored in an `npy` file that contains the following fields:

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

Each cube can be selected from the following drop-down menu:
$(@bind data_path_LDSS3 pl.Select(sort(glob("data/reduced/LDSS3/*/LC*.npy"))))
"""

# ╔═╡ f2bb15ee-2180-4e0f-b71d-7f9cdc2178ef
md"""
#### Common wavelengths
"""

# ╔═╡ 7185f603-3e57-42db-9665-c215094776ad
md"""
Next we extract the common wavelength grid, along with the target and comparison star flux, and compute the resulting spectrum and light curves:
"""

# ╔═╡ 2724283b-78c1-47b3-ac17-8c92e9e677d9
names_LDSS3 = OrderedDict(
	"WASP50" => "aperture_324_803",
	"c06" => "aperture_153_1117",
	"c15" => "aperture_830_689",
	"c21" => "aperture_28_1189",
)

# ╔═╡ 63f2a355-1f29-4a4c-8f6b-efb2007a8e11
md"""
!!! note
	We convert from Float32 to Float64 to match the format of the IMACS data
"""

# ╔═╡ fc9ad91e-ab8c-45ad-b8c6-065ac0cb01fe


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
LC_IMACS = load_pickle(data_path_IMACS)

# ╔═╡ 44808b97-df11-4aff-9e97-f97987fe9939
LC_LDSS3 = load_npz(data_path_LDSS3, allow_pickle=true)

# ╔═╡ 639f666b-09fc-488b-982b-01a523278cae
f_LDSS3 = let
	fluxes = LC_LDSS3["cubes"]["raw_counts"]
	Dict(k => convert(Matrix{Float64}, v) for (k, v) ∈ fluxes)
end

# ╔═╡ 2e503024-65cb-483d-aac7-364f22823bdc
# Non-zero wavelength ranges for each star
idx_lims = [
	maximum(findfirst.(!=(0.0), eachrow(f))):minimum(findlast.(!=(0.0), eachrow(f)))
	for (_, f) in f_LDSS3
]

# ╔═╡ 993e1e5a-9d78-4042-9aca-bb753af7f647
# Use the intersection for common wavelength range
common_wav_idxs_LDSS3 = ∩(idx_lims...)

# ╔═╡ 7c6a3f88-cdc5-4b56-9bbf-2e9a2fa8ae26
f_target_LDSS3 = f_LDSS3[names_LDSS3["WASP50"]][:, common_wav_idxs_LDSS3]

# ╔═╡ bb2c6085-5cb4-4efc-a322-27fda09fb904
f_comps_LDSS3 = cat(
	[
		f_LDSS3[names_LDSS3[comp_name]][:, common_wav_idxs_LDSS3]
		for comp_name in ["c06", "c15", "c21"]
	]...,
	dims=3
)

# ╔═╡ 31bdc830-dfc9-445a-ba7e-76be7627561a
wav_LDSS3 = LC_LDSS3["spectral"]["wavelength"][common_wav_idxs_LDSS3]

# ╔═╡ e774a20f-2d58-486a-ab71-6bde678b26f8
md"""
## Stellar spectra 🌟

We next take a look at the extracted spectra from each instrument.
"""

# ╔═╡ 2bc50f2d-64b4-45e8-83e6-9ae2a948d76a
md"""
### IMACS
"""

# ╔═╡ 0d18676d-5401-44ab-8a95-c45fa7864115
md"""
### LDSS3
"""

# ╔═╡ 6fd88483-d005-4186-8dd2-82cea767ce90
med_std(A; dims=1) = (median(A, dims=dims), std(A, dims=dims)) .|> vec

# ╔═╡ 1f8f5bd0-20c8-4a52-9dac-4ceba18fcc06
function spec_plot(ax, wav, A; color=:blue, label="")
	μ, σ = med_std(A)
	band(ax, wav, μ .- σ, μ .+ σ, color=(color, 0.25))
	lines!(ax, wav, μ, color=color, label=label)
end

# ╔═╡ e3468c61-782b-4f55-a4a1-9d1883655d11
md"""
## White light curves ⚪

Next, we will extract the white light curves for each instrument. We divide the target WLC by each comparison star to i)minimize common systematics (e.g., air mass, local refractive atmospheric effects), and ii) make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ 731ed0be-a679-4738-a477-56b0ed44338c
md"""
### IMACS

First we compute the median normalized, comparison star divided WLCs stored in the pickle file:
"""

# ╔═╡ 18d58341-0173-4eb1-9f01-cfa893088613
begin
	comp_names_IMACS = LC_IMACS["cNames"]
	sorted_cName_idxs = sortperm(comp_names_IMACS)
	
	f_div_WLC_IMACS = LC_IMACS["oLC"] ./ LC_IMACS["cLC"][:, sorted_cName_idxs]
	f_div_WLC_norm_IMACS = f_div_WLC_IMACS ./ median(f_div_WLC_IMACS, dims=1)
end

# ╔═╡ 941cd721-07d8-4a8f-9d75-42854e6e8edb
md"""
!!! warning
	In general, the comparison stars names (`cNames`) are not stored in alphanumeric order by default. For convenience, we ensure this sorting with `sortperm`, so that the first column corresponds to the first comparison star, the second to the second comparison star, etc.

	TODO: Identify cNames order for LDSS3 data

We next plot these light curves and identified outliers below:
"""

# ╔═╡ ab058d99-ce5f-4ed3-97bd-a62d2f258773
@bind window_width pl.Slider(3:2:21, show_value=true)

# ╔═╡ 4b763b58-862e-4c88-a7c9-fe0b1271c0b4
use_comps_IMACS = ["c06"]

# ╔═╡ 9cf2642e-d436-4078-b4a3-e2a519b6d651
md"""
### LDSS3
"""

# ╔═╡ 08eafc08-b0fb-4c99-9d4e-7fa5085d386c
md"""
Using the extracted spectra from earlier, we can integrate the flux along the binned spectral direction **<LINK TO BINNED SECTION>** to build the white-light curves for each star. We store this in the `ntimes` ``\times`` `ncomps` array, `f_div_wlc`:
"""

# ╔═╡ 4b3b0a83-daaf-4339-86cc-24e0ba95ea85
md"""
### Helper functions
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

# ╔═╡ 109ab11b-cbb3-4c02-9b7a-5d6047883364
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ df46d106-f186-4900-9d3f-b711bc803707
pl.with_terminal() do
	use_comps = use_comps_IMACS
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names_IMACS))
	_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm_IMACS[:, use_comps_idxs], window_width)
	# Because python
	println(bad_idxs .- 1)
	println(use_comps_idxs .- 1)
end

# ╔═╡ e98dee2e-a369-448e-bfe4-8fea0f318fa8
md"""
## Binned light curves 🌈
"""

# ╔═╡ 891add6b-c1c1-4e7c-8a8d-de2421bd6f2d
md"""
### IMACS

```
	
```
"""

# ╔═╡ e844ad0c-c3ce-40cc-840f-7fe6ec454fed
md"""
### LDSS3

Similarly to above, we next integrate the target and comparison star flux over the same wavelength scheme as our IMACS data to build the binned light curves below.
"""

# ╔═╡ 88dd7d0b-c133-4fd2-b942-520b7e3d0265
md"""
#### Define wavelength bins

We define the lower and upper bound for each of the `nbins` wavelength bins that the timeseries data will be integrated over, where `nbins` is the total number of wavelength bins used:
"""

# ╔═╡ cfe0f1df-fc45-409f-a4ee-34b48a99c90c
wbins = readdlm("data/reduced/LDSS3/w50_bins_LDSS3.dat", comments=true)

# ╔═╡ a838a0ec-ee9c-4984-a62c-68810ec731de
binned_wav_idxs_LDSS3 = [
	findfirst(==(wbin[1]), wav_LDSS3):findfirst(==(wbin[2]), wav_LDSS3)
	for wbin in eachrow(wbins)
]

# ╔═╡ ee94bd7d-380a-44c4-a8dd-9fd632a96158
binned_wav_idxs_LDSS3_range = binned_wav_idxs_LDSS3[begin][begin]:binned_wav_idxs_LDSS3[end][end]

# ╔═╡ c51edf69-b43b-4119-b2ad-f28f6de48c92
begin
	f_target_wlc_LDSS3 = sum(f_target_LDSS3[:, binned_wav_idxs_LDSS3_range], dims=2)
	f_comps_wlc_LDSS3 = sum(
		f_comps_LDSS3[:, binned_wav_idxs_LDSS3_range, :], dims=2
	) |> x -> dropdims(x, dims=2)
	f_div_wlc_LDSS3 = f_target_wlc_LDSS3 ./ f_comps_wlc_LDSS3
	f_div_WLC_norm_LDSS3 = f_div_wlc_LDSS3 ./ median(f_div_wlc_LDSS3, dims=1)
end

# ╔═╡ 390cb536-1ca1-410d-a751-c361b7349ea2
_, use_idxs_LDSS3, bad_idxs_LDSS3 = filt_idxs(f_div_WLC_norm_LDSS3, window_width);

# ╔═╡ 4891ae49-0311-415b-bbea-c703a81aaf63
pl.with_terminal() do
	println(bad_idxs_LDSS3 .- 1) # For Python
end

# ╔═╡ b9367db3-a608-425e-9fc9-92b4b8879e46
bad_idxs_LDSS3 .- 1 == [0, 2, 3, 4, 5, 7, 9, 10, 12, 14, 15, 18, 21, 22, 23, 24, 25, 26, 28, 31, 32, 34, 37, 38, 39, 40, 41, 42, 43, 44, 46, 48, 51, 52, 53, 54, 55, 63, 75, 76, 78, 79, 82, 196, 201, 208, 209, 210, 213, 215, 217, 233, 235, 247, 254, 264, 270, 279, 283, 289, 304]

# ╔═╡ f7feb44e-a363-4f8d-bf62-d3541533e4da
md"""
#### Compute binned divided LCs
"""

# ╔═╡ eb4f7a92-a9e7-4bf5-8b1c-bca928fcedde
md"""
We will next compute `oLCw` and `cLCw`, where `oLCw` is an `ntimes` ``\times`` `nbins` matrix that holds the binned target flux, where `ntimes` is the number of timeseries points ``N``. Similarly `cLCw` is an `ntimes` ``\times`` `nbins` ``\times`` `ncomps` matrix that holds the comparison star flux, where `ncomps` is the number of comparison stars:
"""

# ╔═╡ 2b20e471-16e5-4b54-abc5-4308af4e60b6
"""
	sum_flux(A::AbstractMatrix, idx_range::AbstractRange, [dims=2])

For columnar data `A`, returns the sum of each row (for default `dims=2`) along the range of columns specified by `idx_range`.
"""
sum_flux(A, idx_range, dims=2) = sum(view(A, :, idx_range), dims=dims)[:, 1]

# ╔═╡ a6805e63-bbf4-48bc-8e15-be24da9b0348
begin
	nbins = length(binned_wav_idxs_LDSS3)
	ntimes, ncomps = size(f_div_wlc_LDSS3)
	oLCw = Matrix{Float64}(undef, ntimes, nbins)
	cLCw = Array{Float64, 3}(undef, ntimes, nbins, ncomps)
	
	for (j, wav_idxs) in enumerate(binned_wav_idxs_LDSS3)
		oLCw[:, j] .= sum_flux(f_target_LDSS3, wav_idxs)
		for k in 1:ncomps
			cLCw[:, j, k] .= sum_flux(f_comps_LDSS3[:, :, k], wav_idxs)
		end
	end
end

# ╔═╡ f44393ab-2e49-4817-9870-890c9cd556ce
md"""
With `oLCw` and `cLCw` now computed, we next compute `f_norm_w`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\times`` `nbins` ``\times`` `ncomps`, where for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve:
"""

# ╔═╡ a73deaa6-5f45-4920-9796-6aa48e80c3de
md"""
We plot these below for each comparison star division. Move the slider to view the plot for the corresponding comparison star:

comp star $(@bind comp_idx pl.Slider(1:ncomps, show_value=true))

!!! note "Future"
	Ability to interact with sliders completely in the browser coming soon!
"""

# ╔═╡ 90fec9cf-37b0-4bcf-a6df-85a4cdfc511b
begin	
	f_norm_w = Array{Float64}(undef, ntimes, nbins, ncomps)
	for c_i in 1:ncomps
		f_w = oLCw ./ cLCw[:, :, c_i]
		f_norm_w[:, :, c_i] .= f_w ./ median(f_w, dims=1)
	end
	
	offs = reshape(range(0, 0.3, length=nbins), 1, :) # Arbitrary offsets for clarity
	f_norm_w .+= offs
	baselines = ones(size(f_norm_w[:, :, comp_idx])) .+ offs # Reference baselines
end;

# ╔═╡ fbc57d8b-3b1b-44d1-bd7d-0e9749026d4c
let
	fig = Figure(resolution = (400, 700))
	ax = fig[1, 1] = Axis(fig, title = "target / comp $comp_idx")
	
	cmap = reverse(to_colormap(:Spectral_4, nbins))
	
	for (c, f, b) in zip(cmap, eachcol(f_norm_w[:, :, comp_idx]), eachcol(baselines))
		scatter!(f, color=c, strokewidth=0, markersize=5)
		lines!(b, color=c)
	end
		
	ylims!(ax, 0.95, 1.34)
	ax.xlabel = "Index"
	ax.ylabel = "Relative flux + offset"
	
	fig |> pl.as_svg
end

# ╔═╡ 12bc210a-03fc-4a85-b281-f129b1877403
md"""
## GPTS inputs 🔩
"""

# ╔═╡ 45268530-f48a-4fbe-89f7-4a7ea972d2b2
md"""
### IMACS

Automatically computed from pickle file
"""

# ╔═╡ 1306cf62-857e-4cb3-a841-d87d2d5a995f
md"""
### LDSS3

Done manually
"""

# ╔═╡ dff24ff0-1d5f-4561-b23f-9385c9e73a0c
md"""
#### External parameters
"""

# ╔═╡ aad465ac-dca8-41a7-a232-c546319dd338
target_name_LDSS3 = names_LDSS3["WASP50"]

# ╔═╡ 9b310d6f-b320-4406-a97c-a32620257995
Times_LDSS3, Airmass_LDSS3 = getindex.(
	Ref(LC_LDSS3["temporal"].columns),
	["bjd", "airmass"]
)

# ╔═╡ 0caa6353-ee39-4e30-b7b0-c5302f87356c
specshifts_LDSS3 = load_npz(
	"$(dirname(data_path_LDSS3))/specshifts.npy", allow_pickle=true
);

# ╔═╡ ec83cd01-ad3b-4c46-a763-828b3f1e0c70
Delta_Wav_LDSS3 = specshifts_LDSS3["shift"][target_name_LDSS3] |> 
		sort |> values |> collect |> x -> convert(Vector{Float64}, x)

# ╔═╡ a77a890a-ddaa-41f4-bdd1-1fef4e56b0f1
fmt_float(x) = @sprintf "%.10f" x

# ╔═╡ 29c61b96-cf27-432f-b5cd-a6192732a8f6
md"""
#### WLCs (magnitude space)
"""

# ╔═╡ abf0021f-6a1e-4a6c-8135-88a0424c3df9
comps_LDSS3 = let
	mag = -2.51 * log10.(f_comps_wlc_LDSS3[use_idxs_LDSS3, :])
	mag .- median(mag, dims=1)
end

# ╔═╡ 9f401644-bf05-4221-8f69-5742785a79b2
md"""
#### Binned LCs (magnitude space)
"""

# ╔═╡ 66a04e1c-8fe1-4aac-95b4-7c126f7d0648
md"""
### Helper functions
"""

# ╔═╡ 62f0c610-862c-466b-b6ba-387ebc11928c
median_eparam(param, cube, obj_name, wav_idxs) = cube[param][obj_name][:, wav_idxs] |>
	x -> median(x, dims=2) |> vec |> x -> convert(Vector{Float64}, x)

# ╔═╡ 0dbc118e-943a-479e-9151-495455d9eed7
FWHM_LDSS3, Trace_Center_LDSS3, Sky_Flux_LDSS3 = median_eparam.(
	["width", "peak", "sky"],
	Ref(LC_LDSS3["cubes"]),
	Ref(target_name_LDSS3),
	Ref(common_wav_idxs_LDSS3[binned_wav_idxs_LDSS3_range])
)

# ╔═╡ 673e676d-a7f2-4b65-8642-8bd508a61bf0
eparams_LDSS3 = DataFrame(
	(Times=Times_LDSS3, Airmass=Airmass_LDSS3, Delta_Wav=Delta_Wav_LDSS3,
	 FWHM=FWHM_LDSS3, Sky_Flux=Sky_Flux_LDSS3, Trace_Center=Trace_Center_LDSS3)
)[use_idxs_LDSS3, :] |> df -> mapcols(col -> fmt_float.(col), df)

# ╔═╡ 655247ce-4c50-4a6c-8186-7dad1233cd3b
CSV.write("$(dirname(data_path_LDSS3))/eparams.dat", eparams_LDSS3, delim=",    ")

# ╔═╡ d127c98d-81c0-4133-955f-ae98c07e1152
function f_to_med_mag(f)
	mag = -2.51 * log10.(f)
	return mag .- median(mag)
end

# ╔═╡ 0e045efd-e9d7-4f03-8128-953c7ad13aba
lc_LDSS3 = let
	med_mag = f_to_med_mag(f_target_wlc_LDSS3 |> vec)
	hcat(Times_LDSS3, med_mag, zeros(length(Times_LDSS3)))[use_idxs_LDSS3, :]
end

# ╔═╡ e282e1da-0389-42f3-aeef-e732f5f1df95
let
	savepath = "$(dirname(data_path_LDSS3))/white-light"
	rm(savepath, force=true, recursive=true)
	mkdir(savepath)
	writedlm("$(savepath)/lc.dat", lc_LDSS3, ",    ")
	writedlm("$(savepath)/comps.dat", comps_LDSS3, ",    ")
end

# ╔═╡ 8576c9d9-29ae-4351-985a-c191a944a4bc
target_binned_mags = mapslices(f_to_med_mag, oLCw, dims=1)[use_idxs_LDSS3, :]

# ╔═╡ 0a7b6bc2-ed75-49d5-b5ad-ff13f86eae2b
comp_binned_mags = mapslices(f_to_med_mag, cLCw, dims=1)[use_idxs_LDSS3, :, :]

# ╔═╡ cb5ced41-cc6b-48f4-a532-0a0f0fe31f8c
let
	savepath = "$(dirname(data_path_LDSS3))/wavelength"
	rm(savepath, force=true, recursive=true)
	for i in 1:nbins
		save_path_w = "$(savepath)/wbin$(i-1)"
		mkpath("$(save_path_w)")
		writedlm("$(save_path_w)/lc.dat", target_binned_mags[:, i], ",    ")
		writedlm("$(save_path_w)/comps.dat", comp_binned_mags[:, i, :], ",    ")
	end
end

# ╔═╡ 5db4a2f2-1c0d-495a-8688-40fc9e0ccd02
md"""
## Plot configs
"""

# ╔═╡ 2b50e77f-0606-4e13-9d8a-c6eb3645d23c
const FIG_TALL = (900, 1_200)

# ╔═╡ 5e97c57e-a896-4478-a277-e3da1444f8de
const FIG_WIDE = (1_350, 800)

# ╔═╡ 202667da-d22b-405d-9935-4726c7d41a0b
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

# ╔═╡ 589239fb-319c-40c2-af16-19025e7b28a2
let
	fig = Figure()
	
	wav = LC_IMACS["spectra"]["wavelengths"]
	spec_plot(fig[1, 1], wav, LC_IMACS["spectra"]["WASP50"];
		color = (COLORS[1]),
	)
	
	fig |> pl.as_svg
end

# ╔═╡ 2fd7bc68-a6ec-4ccb-ad73-07b031ffef5a
let
	fig = Figure(resolution=FIG_WIDE)

	fluxes = [f_target_LDSS3, [f_comps_LDSS3[:, :, i] for i in 1:3]...]
	labels = keys(names_LDSS3) |> collect

	k = 1
	for j in 1:2, i in 1:2
		spec_plot(fig[i, j], wav_LDSS3, fluxes[k];
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

	fig[1:2, 0] = Label(fig, "Counts", rotation=π/2)
	fig[end+1, 2:3] = Label(fig, "Index")
   #textsize = 30, font = noto_sans_bold, color = (:black, 0.25))

	fig |> pl.as_svg
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
	
	ncomps = length(comp_names_IMACS)
	use_comps = use_comps_IMACS
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names_IMACS))
	
	axs = [Axis(fig[i, j]) for i ∈ 1:4, j ∈ 1:2]
	axs = reshape(copy(fig.content), 4, 2)
	
	plot_div_WLCS!(
		axs, f_div_WLC_norm_IMACS, window_width, comp_names_IMACS, use_comps_idxs
	)
	
	linkaxes!(axs...)
	hidexdecorations!.(axs[begin:end-1, :], grid=false)
	hideydecorations!.(axs[:, begin+1:end], grid=false)
	ylims!(axs[end], 0.97, 1.02)
	
	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2)
	fig[end+1, 1:end] = Label(fig, "Index")
	
	fig |> pl.as_svg
end

# ╔═╡ b1bd886f-c7dd-4167-b9b8-084b73f7ee9c
let
	fig = Figure(resolution=FIG_TALL)
	
	comp_names_LDSS3 = collect(keys(names_LDSS3))[begin+1:end]
	ncomps = length(comp_names_LDSS3)
	use_comps = comp_names_LDSS3 # use all comps
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names_LDSS3))
	
	axs = [Axis(fig[i, 1]) for i in 1:ncomps]
	axs = reshape(copy(fig.content), ncomps, 1)
	
	plot_div_WLCS!(
		axs, f_div_WLC_norm_LDSS3, window_width, comp_names_LDSS3, use_comps_idxs
	)
	
	hidexdecorations!.(axs[begin:end-1], grid=false)
	linkaxes!(axs...)
	ylims!(axs[end], 0.97, 1.02)
	
	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2)
	fig[end+1, 1:end] = Label(fig, "Index")
	
	fig |> pl.as_svg
end

# ╔═╡ eeb3da97-72d5-4317-acb9-d28637a06d67
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╟─ee24f7df-c4db-4065-afe9-10be80cbcd6b
# ╟─9d180c21-e634-4a1e-8430-bdd089262f66
# ╟─454cab9c-7ea0-4b0a-9622-d29d8ba0395b
# ╠═bd2cdf33-0c41-4948-82ab-9a28929f72b3
# ╟─66052b03-35a0-4877-abef-f525766fd332
# ╟─470c738f-c2c8-4f56-b936-e08c43c161b1
# ╠═44808b97-df11-4aff-9e97-f97987fe9939
# ╟─f2bb15ee-2180-4e0f-b71d-7f9cdc2178ef
# ╟─7185f603-3e57-42db-9665-c215094776ad
# ╠═2724283b-78c1-47b3-ac17-8c92e9e677d9
# ╠═639f666b-09fc-488b-982b-01a523278cae
# ╟─63f2a355-1f29-4a4c-8f6b-efb2007a8e11
# ╠═2e503024-65cb-483d-aac7-364f22823bdc
# ╠═993e1e5a-9d78-4042-9aca-bb753af7f647
# ╠═31bdc830-dfc9-445a-ba7e-76be7627561a
# ╠═fc9ad91e-ab8c-45ad-b8c6-065ac0cb01fe
# ╠═7c6a3f88-cdc5-4b56-9bbf-2e9a2fa8ae26
# ╟─7bfc971c-8737-49ad-adec-ac57d176f10e
# ╠═bb2c6085-5cb4-4efc-a322-27fda09fb904
# ╟─1c3e8cb3-2eff-47c2-8c17-01d0599556b8
# ╠═3653ee36-35a6-4e0a-8d46-4f8389381d45
# ╟─e774a20f-2d58-486a-ab71-6bde678b26f8
# ╟─2bc50f2d-64b4-45e8-83e6-9ae2a948d76a
# ╠═589239fb-319c-40c2-af16-19025e7b28a2
# ╟─0d18676d-5401-44ab-8a95-c45fa7864115
# ╠═2fd7bc68-a6ec-4ccb-ad73-07b031ffef5a
# ╠═1f8f5bd0-20c8-4a52-9dac-4ceba18fcc06
# ╠═6fd88483-d005-4186-8dd2-82cea767ce90
# ╟─e3468c61-782b-4f55-a4a1-9d1883655d11
# ╟─731ed0be-a679-4738-a477-56b0ed44338c
# ╠═18d58341-0173-4eb1-9f01-cfa893088613
# ╟─941cd721-07d8-4a8f-9d75-42854e6e8edb
# ╠═ab058d99-ce5f-4ed3-97bd-a62d2f258773
# ╠═13523326-a5f2-480d-9961-d23cd51192b8
# ╠═4b763b58-862e-4c88-a7c9-fe0b1271c0b4
# ╠═df46d106-f186-4900-9d3f-b711bc803707
# ╟─9cf2642e-d436-4078-b4a3-e2a519b6d651
# ╟─08eafc08-b0fb-4c99-9d4e-7fa5085d386c
# ╠═c51edf69-b43b-4119-b2ad-f28f6de48c92
# ╠═ee94bd7d-380a-44c4-a8dd-9fd632a96158
# ╠═b1bd886f-c7dd-4167-b9b8-084b73f7ee9c
# ╠═390cb536-1ca1-410d-a751-c361b7349ea2
# ╠═4891ae49-0311-415b-bbea-c703a81aaf63
# ╠═b9367db3-a608-425e-9fc9-92b4b8879e46
# ╟─4b3b0a83-daaf-4339-86cc-24e0ba95ea85
# ╠═ccabf5d2-5739-4284-a972-23c02a263a5c
# ╠═a4517d69-76e6-462a-9449-b31d80e34a8f
# ╠═dc044a72-4706-49e2-94a8-c828a6bf7de0
# ╠═109ab11b-cbb3-4c02-9b7a-5d6047883364
# ╟─e98dee2e-a369-448e-bfe4-8fea0f318fa8
# ╟─891add6b-c1c1-4e7c-8a8d-de2421bd6f2d
# ╟─e844ad0c-c3ce-40cc-840f-7fe6ec454fed
# ╟─88dd7d0b-c133-4fd2-b942-520b7e3d0265
# ╠═cfe0f1df-fc45-409f-a4ee-34b48a99c90c
# ╠═a838a0ec-ee9c-4984-a62c-68810ec731de
# ╟─f7feb44e-a363-4f8d-bf62-d3541533e4da
# ╟─eb4f7a92-a9e7-4bf5-8b1c-bca928fcedde
# ╠═2b20e471-16e5-4b54-abc5-4308af4e60b6
# ╠═a6805e63-bbf4-48bc-8e15-be24da9b0348
# ╟─f44393ab-2e49-4817-9870-890c9cd556ce
# ╠═90fec9cf-37b0-4bcf-a6df-85a4cdfc511b
# ╟─a73deaa6-5f45-4920-9796-6aa48e80c3de
# ╠═fbc57d8b-3b1b-44d1-bd7d-0e9749026d4c
# ╟─12bc210a-03fc-4a85-b281-f129b1877403
# ╟─45268530-f48a-4fbe-89f7-4a7ea972d2b2
# ╟─1306cf62-857e-4cb3-a841-d87d2d5a995f
# ╟─dff24ff0-1d5f-4561-b23f-9385c9e73a0c
# ╠═aad465ac-dca8-41a7-a232-c546319dd338
# ╠═9b310d6f-b320-4406-a97c-a32620257995
# ╠═0dbc118e-943a-479e-9151-495455d9eed7
# ╠═0caa6353-ee39-4e30-b7b0-c5302f87356c
# ╠═ec83cd01-ad3b-4c46-a763-828b3f1e0c70
# ╠═673e676d-a7f2-4b65-8642-8bd508a61bf0
# ╠═a77a890a-ddaa-41f4-bdd1-1fef4e56b0f1
# ╠═442946ba-9561-4223-8223-ba513ac7bc45
# ╠═655247ce-4c50-4a6c-8186-7dad1233cd3b
# ╟─29c61b96-cf27-432f-b5cd-a6192732a8f6
# ╠═0e045efd-e9d7-4f03-8128-953c7ad13aba
# ╠═abf0021f-6a1e-4a6c-8135-88a0424c3df9
# ╠═e282e1da-0389-42f3-aeef-e732f5f1df95
# ╟─9f401644-bf05-4221-8f69-5742785a79b2
# ╠═8576c9d9-29ae-4351-985a-c191a944a4bc
# ╠═0a7b6bc2-ed75-49d5-b5ad-ff13f86eae2b
# ╠═cb5ced41-cc6b-48f4-a532-0a0f0fe31f8c
# ╟─66a04e1c-8fe1-4aac-95b4-7c126f7d0648
# ╠═62f0c610-862c-466b-b6ba-387ebc11928c
# ╠═d127c98d-81c0-4133-955f-ae98c07e1152
# ╟─5db4a2f2-1c0d-495a-8688-40fc9e0ccd02
# ╠═2b50e77f-0606-4e13-9d8a-c6eb3645d23c
# ╠═5e97c57e-a896-4478-a277-e3da1444f8de
# ╠═202667da-d22b-405d-9935-4726c7d41a0b
# ╟─eeb3da97-72d5-4317-acb9-d28637a06d67
# ╠═b1b0690a-a1eb-11eb-1590-396d92c80c23
