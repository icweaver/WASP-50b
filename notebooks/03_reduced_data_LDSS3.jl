### A Pluto.jl notebook ###
# v0.17.5

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
	using AlgebraOfGraphics
	using CairoMakie
	using CSV, DataFrames
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using Statistics
	using OrderedCollections
	using Printf
end

# ╔═╡ bda7227d-a952-4380-ade4-7cf784a1e5cd
using Dates

# ╔═╡ 26f18ff6-7baa-4905-b9d1-52cfa9396dfc
begin
	ENV["PYTHON"] = ""
	Pkg.build("PyCall")
	using PyCall
	const Conda = PyCall.Conda
	Conda.add("lightkurve", :WASP50b)
end

# ╔═╡ 34ef4580-bb95-11eb-34c1-25893217f422
md"""
# Reduced Data -- LDSS3

In this notebook we will examine the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from LDSS3.

$(TableOfContents(depth=4))
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
"""

# ╔═╡ 48713a10-6ba7-4de5-a147-85a9cb136dd8
@bind FPATH Select(glob("data/reduced_data/LDSS3/ut150927_*/LCs_*.npy"))

# ╔═╡ 99a310b2-b802-4b29-9354-796b2dec2820
const FIG_PATH = "figures/reduced_data"

# ╔═╡ 4b8d9fc4-3d00-4311-89c4-56478122b2b4
julian2datetime.((2457292.705208333,2457292.7103152038))

# ╔═╡ 698ae5b0-7cd3-4055-a13e-e9aa3704ca12
md"""
We next define the mapping between the standard aperture names defined in `stellar` and the target/comparison star names used:
"""

# ╔═╡ b4f4e65b-bd50-4ee2-945f-7c130db21fdf
md"""
!!! note
	These stars were cross-matched with VizieR to verify the aperture name assignments

	**Add to paper:** "This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: 10.26093/cds/vizier). The original description of the VizieR service was published in A&AS 143, 23"
"""

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
## $(@bind plot_stell_spec CheckBox()) Stellar spectra ⭐

With the flux extracted for each object, we now turn to analyzing the resulting stellar spectra, selected from the wavelength bins scheme below:
"""

# ╔═╡ 9697e26b-b6d9-413b-869f-47bc2ab99919
@bind FPATH_WBINS let
	dirpath = "data/reduced_data/wbins"
	Select(["$(dirpath)/w50_bins_LDSS3.dat", "$(dirpath)/w50_bins_species.dat"])
end

# ╔═╡ db933939-b8df-48b2-b53e-1dbd7ec1b07c
wbins = readdlm(FPATH_WBINS, comments=true)

# ╔═╡ 25903146-2f23-4e54-bd7d-aed1025f2fa5
fname_suff = let
	suff = "LDSS3_" * basename(dirname(FPATH))
	occursin("species", FPATH_WBINS) ? (suff *= "_species") : suff
end

# ╔═╡ 979895ea-9b13-494c-be08-2496562ccf07
Na_bins = range(5800.0, stop=5986.0, length=6)

# ╔═╡ ce86c8d1-1900-413f-90be-c9f4eb386a5a
K_bins = range(7655.0, stop=7709.0, length=6)

# ╔═╡ bf61788d-ca13-457e-8026-62bc2460d4d7
make_bins(x) = hcat(x[begin:end-1], x[begin+1:end])

# ╔═╡ 48f31cef-ee18-47f9-a71d-4a5b2c6efe71
#species = (Na=5892.9, K=7682.0, Na_8200=8189.0)

# ╔═╡ a6088ea2-904f-4909-b1be-9470e7ec2010
med_std(A; dims=1) = (median(A, dims=dims), std(A, dims=dims)) .|> vec

# ╔═╡ 7d68ad39-3e39-48fa-939a-e56c6659d2b3
function spec_plot!(ax, wav, A; color=:blue, norm=1.0, label="")
	μ, σ = med_std(A) ./ norm
	band!(ax, wav, μ .- σ, μ .+ σ, color=(color, 0.25))
	lines!(ax, wav, μ;
		color = color,
		cycle = Cycle(:linestyle),
		label = label,
	)
end

# ╔═╡ bd937d51-17e9-4de3-a5d0-4c436d413940
md"""
## $(@bind plot_lcs CheckBox()) White-light curves 🌅

Next, we will extract the integrated white-light curves from these spectra. We integrate over the same wavelength bins used in the IMACS analysis:
"""

# ╔═╡ cb805821-5d2e-484c-93a5-10a897d2cfe7
@bind window_width PlutoUI.Slider(3:2:21, default=15, show_value=true)

# ╔═╡ 9372c69a-0aad-4e6e-9ea3-e934fa09b758
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ d5c6d058-17c6-4cf0-97b8-d863b1529161
md"""
!!! note
	We divide the target WLC by each comparison star to minimize common systematics (e.g., air mass, local refractive atmospheric effects), and to make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ 3dc704cc-abee-42b8-8220-62fdda278944
xs = [0.5, 0.4, 0.8, 0.6, 0.55, 0.7, 0.45]

# ╔═╡ 3d4dd810-e317-40da-8521-365c7b760fe4
let
	fig = Figure()
	ax = Axis(fig[1, 1])

	scatter!(ax, 1:7, xs)
	errorbars!(ax, 1:7, xs, std(xs))
	hlines!(ax, median(xs), color=:red)
	
	fig
end

# ╔═╡ 8d91fe42-cd8d-4ad1-bed4-4cee3b778cfd
xs

# ╔═╡ ae33702b-700c-4a4b-83b0-c9f36bf1404f
mapwindow(median, xs, 5)

# ╔═╡ 89256633-3c12-4b94-b245-4fdda44d686c
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
md"""
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
	
# 	comp_names = obj_names.vals[2:4]
# 	f_wlc_targ = f_target_wlc[:]
# 	f_wlc_comps = f_comps_wlc ./ mean(f_comps_wlc, dims=1)
	
# 	for (cName, col) in zip(sort(comp_names), eachcol(f_wlc_comps))
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
md"""
## $(@bind plot_blcs CheckBox()) Binned light curves 🌈

We next integrate the target and comparison star flux within each bin defined above to build the binned light curves. We store these in `oLCw` and `cLCw`, where `oLCw` is an `ntimes` ``\times`` `nbins` matrix that holds the binned target flux, where `ntimes` is the number of timeseries points ``N``. Similarly `cLCw` is an `ntimes` ``\times`` `ncomps` ``\times`` `nbins` matrix that holds the comparison star flux, where `ncomps` is the number of comparison stars:
"""

# ╔═╡ 90b26a58-834a-445b-9242-d7b2e04be614
# For columnar data `A`, returns the sum of each row (for default `dims=2`)
# along the range of columns specified by `idx_range`
sum_flux(A, idx_range, dims=2) = sum(view(A, :, idx_range), dims=dims)[:, 1]

# ╔═╡ 5837dc0e-3537-4a90-bd2c-494e7b3e6bb7
md"""
With `oLCw` and `cLCw` now computed, we next compute `f_norm_w`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\times`` `ncomps` ``\times`` `nbins`, where for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve:
"""

# ╔═╡ af07dc54-eeb5-4fbe-8dd0-289dea97502a
md"""
## `GPTransmissionSpectra` inputs 🔩

Finally, we export the external parameters and light curves (in magnitude space) to be used for detrending in `GPTransmissionSpectra`:
"""

# ╔═╡ 88cc640d-b58d-4cde-b793-6c66e74f6b3a
md"""
#### External parameters
"""

# ╔═╡ 301ff07c-8dd5-403a-bae8-a4c38deeb331
target_name = "aperture_324_803" #"aperture_325_803" # "aperture_324_803"

# ╔═╡ d14ab9de-23b4-4647-a823-9b318bb734e9
median_eparam(param, cube, obj_name, wav_idxs) = cube[param][obj_name][:, wav_idxs] |>
	x -> median(x, dims=2) |> vec |> x -> convert(Vector{Float64}, x)

# ╔═╡ e6ed8da6-189a-4931-aac5-a4e8d0291723
fmt_float(x) = @sprintf "%.10f" x

# ╔═╡ 2aba612a-7599-4a2d-9ff0-2fd398c2a0db
let
    #savepath = template_dir(FPATH)
	#rm(savepath, force=true, recursive=true)
	#mkpath(savepath)
end

# ╔═╡ 0be35b52-caea-4000-8cf8-ab99205bdb97
function template_dir(fpath)
	base_dir = "$(dirname(dirname(fpath)))/LDSS3_template/WASP50"
	date, flat_status = split(basename(dirname(fpath)), '_')
	return "$(base_dir)/w50_$(date[3:end])_LDSS3_$(flat_status)"
end

# ╔═╡ 4af35a95-99b3-4186-a47d-169b9cbf927b
tdir = "$(template_dir(FPATH))"

# ╔═╡ 2341920d-6db8-4f08-b8ed-4907ceab7357
md"""
#### WLCs (magnitude space)
"""

# ╔═╡ 898a390b-49f7-45f4-b1a1-b22922d69a29
let
    #savepath = "$(tdir)/white-light"
	#rm(savepath, force=true, recursive=true)
	#mkpath(savepath)
	#writedlm("$(savepath)/lc.dat", lc, ",    ")
	#writedlm("$(savepath)/comps.dat", comps, ",    ")
end

# ╔═╡ b5ef9d95-1a2a-4828-8e27-89727f2e288b
md"""
#### Binned LCs (magnitude space)
"""

# ╔═╡ 631c4d02-58cc-4c70-947f-22c8e8b11015
let
    #savepath = "$(tdir)/wavelength"
	#rm(savepath, force=true, recursive=true)
	#for i in 1:nbins
	#	save_path_w = "$(savepath)/wbin$(i-1)"
	#	mkpath("$(save_path_w)")
	#	lc_w = hcat(
	#		times[use_idxs], target_binned_mags[:, i], zeros(length(times[use_idxs]))
	#	)
	#	writedlm("$(save_path_w)/lc.dat", lc_w, ",    ")
	#	writedlm("$(save_path_w)/comps.dat", comp_binned_mags[:, :, i], ",    ")
	#end
end

# ╔═╡ 123d0c63-f05a-4a7d-be16-6a3b9abac044
function f_to_med_mag(f)
	mag = -2.51 * log10.(f)
	return mag .- median(mag)
end

# ╔═╡ f788835c-8e81-4afe-805e-4caf2d5e5d5b
md"""
## Notebook setup
"""

# ╔═╡ 36c1aa6d-cde9-4ff0-b55c-13f43e94256d
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 2d34e125-548e-41bb-a530-ba212c0ca17c
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
LC = load_npz(FPATH, allow_pickle=true)

# ╔═╡ d8236985-9a36-4357-ac78-7eb39dd0f080
obj_names = OrderedDict(
	LC["target"] => "WASP50",
	"aperture_153_1117" => "c06",
	"aperture_830_689" => "c15",
	"aperture_28_1189" => "c21",
	"aperture_157_1116" => "c06",
	"aperture_830_693" => "c15",
	"aperture_35_1188" => "c21",
)

# ╔═╡ 83263daf-a902-4414-850b-aa6949752fbb
comp_names = obj_names.vals[2:4]

# ╔═╡ cd10cbf3-22f4-46cd-8345-cec3d141e3ca
use_comps = comp_names 

# ╔═╡ 4f71ba8d-bfa0-4adc-8d82-cd3bca8b6c14
use_comps_idxs = get_idx.(use_comps, Ref(comp_names))

# ╔═╡ 9341c427-642a-4782-925a-6e07b91277a0
fluxes = Dict(obj_names[k] => v for (k, v) ∈ LC["cubes"]["raw_counts"])

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

# ╔═╡ ded63b4b-61b6-41b6-98d4-d13166bce76a
@with_terminal begin
	println(bad_idxs .- 1) # For Python
end

# ╔═╡ 3e7b0a0b-1ee9-4436-935e-c4ced50620ba
size(f_div_WLC_norm, 1) - length(bad_idxs)

# ╔═╡ bc54942e-38ef-4919-a3d4-28d5f4db8487
comps = let
	mag = -2.51 * log10.(f_comps_wlc)
	mag .-= median(mag, dims=1)
	mag[use_idxs, :]
end

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

# ╔═╡ 5110de9a-3721-4043-b8b7-493daacb4137
target_binned_mags = mapslices(f_to_med_mag, oLCw, dims=1)[use_idxs, :]

# ╔═╡ 53f5a645-93e0-499a-bb36-e4ff0153a63c
comp_binned_mags = mapslices(f_to_med_mag, cLCw, dims=1)[use_idxs, :, :]

# ╔═╡ c03cb527-d16d-47aa-ab63-6970f4ff0b1f
times, airmass = let
	vals = LC["temporal"].columns
	vals["bjd"], vals["airmass"]
end

# ╔═╡ 2f8bbd81-2f4b-4ca3-bf01-a27bca3a0f19
t = times

# ╔═╡ 32fb9623-a711-4582-aac9-38f46513d742
round.(julian2datetime.((t[begin], t[end])), Dates.Minute)

# ╔═╡ 2bc4f498-9732-459f-9d8e-c1c549562c63
extrema(airmass)

# ╔═╡ 9cdee207-8911-4d37-a5bd-b920e5a8846b
airmass[[begin, end]]

# ╔═╡ 523aed12-afad-408a-a13f-85ba1b86315d
lines(airmass)

# ╔═╡ 354580e4-0aa9-496f-b024-665025a2eeda
lc = let
	med_mag = f_to_med_mag(f_target_wlc |> vec)
	hcat(times, med_mag, zeros(length(times)))[use_idxs, :]
end

# ╔═╡ 1c3e3620-e4de-401d-afee-66303d35a9e2
LC["temporal"].columns["exptime"] |> unique

# ╔═╡ c65c2298-e3a3-4666-be9d-73ee43d94847
fwhm, trace_center, sky_flux = median_eparam.(
	["width", "peak", "sky"],
	Ref(LC["cubes"]),
	Ref(target_name),
	Ref(common_wav_idxs[binned_wav_idxs_range])
)

# ╔═╡ 079c3915-33af-40db-a544-28453732c372
specshifts = load_npz(
	"$(dirname(FPATH))/specshifts.npy", allow_pickle=true
);

# ╔═╡ e4960d1a-8e33-478a-8100-d1838782938d
delta_wav = specshifts["shift"][target_name] |> 
		sort |> values |> collect |> x -> convert(Vector{Float64}, x)

# ╔═╡ e4388fba-64ef-4588-a1ed-283da2f52196
df_eparams = DataFrame(
	(Times=times, Airmass=airmass, Delta_Wav=delta_wav,
	 FWHM=fwhm, Sky_Flux=sky_flux, Trace_Center=trace_center)
) |> df -> mapcols(col -> fmt_float.(col), df)[use_idxs, :]

# ╔═╡ af7beb47-ecfa-4d08-b239-c27f7e8bdc4c
CSV.write("$(tdir)/eparams.dat", df_eparams, delim=",    ")

# ╔═╡ c911cecd-0747-4cd1-826f-941f2f58091c
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

# ╔═╡ 45418bd3-74a3-4758-9fce-adddbeeec076
if plot_stell_spec let
	fig = Figure(resolution=FIG_WIDE)
	ax = Axis(fig[1, 1];
		xlabel = "Wavelength Å",
		ylabel = "Relative flux",
	)
	
	fluxes = [f_target, [f_comps[:, :, i] for i in 1:3]...]
	labels = obj_names.vals
	#f_norm = 40362.188283796 # median IMACS WASP-50 flux for comparison
	
	for (i, (name, f)) in enumerate(zip(labels, fluxes))
		spec_plot!(ax, wav, f;
			color=COLORS_SERIES[i],
			#norm = f_norm,
			norm = median(f),
			label = name,
		)
	end

	vlines!.(ax, wbins, linewidth=1.0, color=:lightgrey)
	
	xlims!(ax, 4_500, 11_000)
	ylims!(ax, 0, 2.6)
	
	axislegend("Transit 2 (LDSS3)")
	
	savefig(fig, "$(FIG_PATH)/extracted_spectra_$(fname_suff).png")
	
	fig
	end
end

# ╔═╡ 20d12d7b-c666-46c3-8f48-5501641e8df3
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

# ╔═╡ 4b2ed9db-0a17-4e52-a04d-3a6a5bf2c054
if plot_lcs let
	fig = Figure(resolution=FIG_WIDE)
	
	# comp_names = obj_names.vals[2:4]
	# ncomps = length(comp_names)
	# use_comps = comp_names # use all comps
	# use_comps_idxs = get_idx.(use_comps, Ref(comp_names))
	
	axs = [Axis(fig[1, i]) for i in 1:ncomps]
	axs = reshape(copy(fig.content), ncomps, 1)
	
	plot_div_WLCS!(
		axs, f_div_WLC_norm, window_width, comp_names, use_comps_idxs
	)
	
	hideydecorations!.(axs[begin+1:end], grid=false)
	linkaxes!(axs...)
	ylims!(axs[end], 0.97, 1.02)
	
	fig[:, 0] = Label(fig, "Relative flux", rotation=π/2, tellheight=false)
	axs[2].xlabel = "Index"

	Label(fig[0, end], "Transit 2 (LDSS3)";
		tellwidth = false,
		halign = :right,
		font = AlgebraOfGraphics.firasans("Bold"),
	)

	savefig(fig, "$(FIG_PATH)/div_wlcs_$(fname_suff).png")
	
	fig
	end
end

# ╔═╡ 31b7e180-207e-459b-bd6a-a965da0ad70c
function plot_binned_lcs(comp_idx; show_fig=false)
	fig = Figure(resolution=FIG_TALL)
	
	comp_name = obj_names.vals[begin+1:end][comp_idx]
	ax_left = Axis(fig[1, 1], title = "target / $(comp_name)")
	ax_right = Axis(fig[1, 2], title = "residuals")
	ax_label = Axis(fig[1, 3])
	axs = reshape(copy(fig.content), 1, 3)
	linkaxes!(axs...)

	# Median filtered curves for visualization purposes
	f_med, _, f_diff = filt(f_norm_w[:, comp_idx, :], window_width)
		
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

	savefig(fig, "$(FIG_PATH)/div_blcs_$(fname_suff)_$(comp_name).png")

	fig
end

# ╔═╡ 9a6d25a2-6a44-49a7-a9e8-aa3651d67ae0
plot_blcs && plot_binned_lcs.(use_comps_idxs)

# ╔═╡ 0cd62849-c726-47d3-94dd-625e6c058cb1
@with_terminal Conda.list(:WASP50b)

# ╔═╡ b2c61d08-6fcf-4b0c-a21a-c0c5e3205210
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
# ╟─34ef4580-bb95-11eb-34c1-25893217f422
# ╟─a8e3b170-3fc2-4d12-b117-07bd37d27710
# ╟─48713a10-6ba7-4de5-a147-85a9cb136dd8
# ╟─99a310b2-b802-4b29-9354-796b2dec2820
# ╠═01e27cf6-0a1b-4741-828a-ef8bf7037ae9
# ╠═4b8d9fc4-3d00-4311-89c4-56478122b2b4
# ╠═32fb9623-a711-4582-aac9-38f46513d742
# ╟─698ae5b0-7cd3-4055-a13e-e9aa3704ca12
# ╠═d8236985-9a36-4357-ac78-7eb39dd0f080
# ╟─b4f4e65b-bd50-4ee2-945f-7c130db21fdf
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
# ╟─9697e26b-b6d9-413b-869f-47bc2ab99919
# ╟─db933939-b8df-48b2-b53e-1dbd7ec1b07c
# ╟─25903146-2f23-4e54-bd7d-aed1025f2fa5
# ╠═45418bd3-74a3-4758-9fce-adddbeeec076
# ╠═979895ea-9b13-494c-be08-2496562ccf07
# ╠═ce86c8d1-1900-413f-90be-c9f4eb386a5a
# ╠═bf61788d-ca13-457e-8026-62bc2460d4d7
# ╠═48f31cef-ee18-47f9-a71d-4a5b2c6efe71
# ╠═7d68ad39-3e39-48fa-939a-e56c6659d2b3
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
# ╠═3e7b0a0b-1ee9-4436-935e-c4ced50620ba
# ╠═4b2ed9db-0a17-4e52-a04d-3a6a5bf2c054
# ╠═9372c69a-0aad-4e6e-9ea3-e934fa09b758
# ╟─d5c6d058-17c6-4cf0-97b8-d863b1529161
# ╠═3dc704cc-abee-42b8-8220-62fdda278944
# ╠═3d4dd810-e317-40da-8521-365c7b760fe4
# ╠═8d91fe42-cd8d-4ad1-bed4-4cee3b778cfd
# ╠═ae33702b-700c-4a4b-83b0-c9f36bf1404f
# ╠═20d12d7b-c666-46c3-8f48-5501641e8df3
# ╠═470514e7-0f08-44a3-8519-5d704ea6b8d4
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
# ╠═9a6d25a2-6a44-49a7-a9e8-aa3651d67ae0
# ╠═31b7e180-207e-459b-bd6a-a965da0ad70c
# ╟─af07dc54-eeb5-4fbe-8dd0-289dea97502a
# ╟─88cc640d-b58d-4cde-b793-6c66e74f6b3a
# ╠═301ff07c-8dd5-403a-bae8-a4c38deeb331
# ╠═c03cb527-d16d-47aa-ab63-6970f4ff0b1f
# ╠═2f8bbd81-2f4b-4ca3-bf01-a27bca3a0f19
# ╠═bda7227d-a952-4380-ade4-7cf784a1e5cd
# ╠═2bc4f498-9732-459f-9d8e-c1c549562c63
# ╠═9cdee207-8911-4d37-a5bd-b920e5a8846b
# ╠═523aed12-afad-408a-a13f-85ba1b86315d
# ╠═1c3e3620-e4de-401d-afee-66303d35a9e2
# ╠═c65c2298-e3a3-4666-be9d-73ee43d94847
# ╠═d14ab9de-23b4-4647-a823-9b318bb734e9
# ╠═079c3915-33af-40db-a544-28453732c372
# ╠═e4960d1a-8e33-478a-8100-d1838782938d
# ╠═e4388fba-64ef-4588-a1ed-283da2f52196
# ╠═e6ed8da6-189a-4931-aac5-a4e8d0291723
# ╠═2aba612a-7599-4a2d-9ff0-2fd398c2a0db
# ╠═4af35a95-99b3-4186-a47d-169b9cbf927b
# ╠═af7beb47-ecfa-4d08-b239-c27f7e8bdc4c
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
# ╠═0cd62849-c726-47d3-94dd-625e6c058cb1
# ╠═26f18ff6-7baa-4905-b9d1-52cfa9396dfc
# ╠═f883b759-65fc-466e-9c8f-e4f941def935
# ╟─b2c61d08-6fcf-4b0c-a21a-c0c5e3205210
