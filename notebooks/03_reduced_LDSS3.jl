### A Pluto.jl notebook ###
# v0.16.3

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

# ╔═╡ c911cecd-0747-4cd1-826f-941f2f58091c
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	Pkg.instantiate()
	
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using CCDReduction
	using Colors
	using DataFrames
	using DataFramesMeta
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using Measurements
	using Measurements: value, uncertainty
	using NaturalSort
	using OrderedCollections
	using Printf
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
	
	import CairoMakie.Makie.KernelDensity: kde
	
	# Python setup
	ENV["PYTHON"] = "/home/mango/miniconda3/envs/WASP-50b/bin/python"
	Pkg.build("PyCall")
	using PyCall
	
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = parse.(Colorant,
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
$(@bind fpath Select(sort(glob("data/reduced/LDSS3/*/LC*.npy"))))
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
	
	
	"aperture_325_803" => "WASP50",
	"aperture_157_1116" => "c06",
	"aperture_830_693" => "c15",
	"aperture_35_1188" => "c21",
)

# ╔═╡ b4f4e65b-bd50-4ee2-945f-7c130db21fdf
md"""
!!! note
	These stars were cross-matched with VizieR to verify the aperture name assignments

	**Add to paper:** "This research has made use of the VizieR catalogue access tool, CDS, Strasbourg, France (DOI: 10.26093/cds/vizier). The original description of the VizieR service was published in A&AS 143, 23"
"""

# ╔═╡ 0d017105-e3a3-4a0e-9e5f-be4200b8cc2c
CSV.File("data/reduced/LDSS3/catalog.tsv", comment="#")

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
## Stellar spectra ⭐

With the flux extracted for each object, we now turn to analyzing the resulting stellar spectra:
"""

# ╔═╡ 979895ea-9b13-494c-be08-2496562ccf07
Na_bins = range(5800.0, stop=5986.0, length=6)

# ╔═╡ ce86c8d1-1900-413f-90be-c9f4eb386a5a
K_bins = range(7655.0, stop=7709.0, length=6)

# ╔═╡ bf61788d-ca13-457e-8026-62bc2460d4d7
make_bins(x) = hcat(x[begin:end-1], x[begin+1:end])

# ╔═╡ 86cedb14-37d1-494b-b022-181887195719
wbins_species = round.(
	vcat(make_bins(Na_bins), make_bins(K_bins)),
)

# ╔═╡ 80d2f86f-377f-4dac-9b89-89518f001a6a
writedlm("./data/reduced/w50_bins_species.dat", wbins_species, " ")

# ╔═╡ 48f31cef-ee18-47f9-a71d-4a5b2c6efe71
species = (Na=5892.9, K=7682.0, Na_8200=8189.0)

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
## White-light curves 🌅

Next, we will extract the integrated white-light curves from these spectra. We integrate over the same wavelength bins used in the IMACS analysis:
"""

# ╔═╡ 9697e26b-b6d9-413b-869f-47bc2ab99919
wbins = readdlm("data/reduced/LDSS3/w50_bins_LDSS3.dat", comments=true)
#wbins = readdlm("data/reduced/w50_bins_species.dat", comments=true)

# ╔═╡ 80480fc2-861b-4fa9-be2d-fcb57f479a93
wbins

# ╔═╡ cb805821-5d2e-484c-93a5-10a897d2cfe7
@bind window_width Slider(3:2:21, default=15, show_value=true)

# ╔═╡ 9372c69a-0aad-4e6e-9ea3-e934fa09b758
get_idx(needle, haystack) = findfirst(==(needle), haystack)

# ╔═╡ d5c6d058-17c6-4cf0-97b8-d863b1529161
md"""
!!! note
	We divide the target WLC by each comparison star to minimize common systematics (e.g., air mass, local refractive atmospheric effects), and to make the transit shape more apparent. This allows us to select good comparison stars for that particular night and which timeseries points to include in the rest of the analysis.
"""

# ╔═╡ fda9022a-2b95-4d06-8a80-80321a1bf441
x = rand(100, 100)

# ╔═╡ 1d91f3b2-293f-467c-ba40-7dc6497152cd
[
 1 2
 3     4	
]

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

# ╔═╡ 97191c55-f673-46b4-82fd-147e7d150623
md"""
### Raw flux
"""

# ╔═╡ a5e742e5-fcca-40d7-b342-c6112e6899e5
md"""
## Binned light curves 🌈

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
target_name = "aperture_324_803"

# ╔═╡ d14ab9de-23b4-4647-a823-9b318bb734e9
median_eparam(param, cube, obj_name, wav_idxs) = cube[param][obj_name][:, wav_idxs] |>
	x -> median(x, dims=2) |> vec |> x -> convert(Vector{Float64}, x)

# ╔═╡ e6ed8da6-189a-4931-aac5-a4e8d0291723
fmt_float(x) = @sprintf "%.10f" x

# ╔═╡ 0be35b52-caea-4000-8cf8-ab99205bdb97
function template_dir(fpath)
	base_dir = "$(dirname(dirname(fpath)))/LDSS3_template/WASP50"
	date, flat_status = split(basename(dirname(fpath)), '_')
	return "$(base_dir)/w50_$(date[3:end])_LDSS3_$(flat_status)"
end

# ╔═╡ 2aba612a-7599-4a2d-9ff0-2fd398c2a0db
let
	savepath = template_dir(fpath)
	rm(savepath, force=true, recursive=true)
	mkpath(savepath)
end

# ╔═╡ 4af35a95-99b3-4186-a47d-169b9cbf927b
tdir = "$(template_dir(fpath))"

# ╔═╡ 2341920d-6db8-4f08-b8ed-4907ceab7357
md"""
#### WLCs (magnitude space)
"""

# ╔═╡ b5ef9d95-1a2a-4828-8e27-89727f2e288b
md"""
#### Binned LCs (magnitude space)
"""

# ╔═╡ 123d0c63-f05a-4a7d-be16-6a3b9abac044
function f_to_med_mag(f)
	mag = -2.51 * log10.(f)
	return mag .- median(mag)
end

# ╔═╡ 3b6b57e2-46ab-46d9-b334-6264daf583f3
md"""
## Helper functions
"""

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
LC = load_npz(fpath, allow_pickle=true)

# ╔═╡ 3ce1d29a-7219-49be-baba-4125b5a557de
am = LC["temporal"].columns["airmass"]

# ╔═╡ de51df33-1be5-473b-afa9-c5a2fc40493a
lines(am)

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

# ╔═╡ f6ed9c49-bf9c-467d-a602-8d18758a4602
fluxes["WASP50"]

# ╔═╡ 2312f321-67ff-4fec-ba73-ff78c82af8dd
median(fluxes["WASP50"], dims=2)

# ╔═╡ dc62887d-746b-4503-8547-7a6814de66a8
wav = LC["spectral"]["wavelength"][common_wav_idxs]

# ╔═╡ 45418bd3-74a3-4758-9fce-adddbeeec076
let
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="Wavelength Å", ylabel="Relative flux")
	
	fluxes = [f_target, [f_comps[:, :, i] for i in 1:3]...]
	labels = names.vals
	#f_norm = 40362.188283796 # median IMACS WASP-50 flux for comparison
	
	for (i, (name, f)) in enumerate(zip(labels, fluxes))
		spec_plot!(ax, wav, f;
			color=COLORS_SERIES[i],
			#norm = f_norm,
			norm = median(f),
			label = name,
		)
	end
	
	#vlines!.(ax, wbins, linewidth=15.0, color=:lightgrey)
	#vlines!.(ax, values(species), color=:red)
	
	#vlines!(ax, Na_bins, color=:green)
	#vlines!(ax, K_bins, color=:green)
	#vlines!.(ax, wbins_species, color=:green)
	
	#xlims!(ax, 5_700, 6_400)
	#xlims!(ax, 7600, 7800)
	
	ylims!(ax, 0, 2.6)
	
	axislegend()
	
	save(
		"../../ACCESS_WASP-50b/figures/reduced/extracted_spectra_ut150927_LDSS3.pdf",
		fig
	)
	
	fig #|> as_svg
end

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

# ╔═╡ 4b2ed9db-0a17-4e52-a04d-3a6a5bf2c054
let
	fig = Figure()#resolution=FIG_TALL)
	
	comp_names = names.vals[2:4]
	ncomps = length(comp_names)
	use_comps = comp_names # use all comps
	use_comps_idxs = get_idx.(use_comps, Ref(comp_names))
	
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
	
	save("../../ACCESS_WASP-50b/figures/reduced/div_wlcs_ut150927_LDSS3.pdf", fig)
	
	fig #|> as_svg
end

# ╔═╡ 470514e7-0f08-44a3-8519-5d704ea6b8d4
_, use_idxs, bad_idxs = filt_idxs(f_div_WLC_norm, window_width);

# ╔═╡ ded63b4b-61b6-41b6-98d4-d13166bce76a
with_terminal() do
	println(bad_idxs .- 1) # For Python
end

# ╔═╡ 56efa000-a943-4256-94f4-f0ae4764f634
let
	fig = Figure()
	ax = Axis(fig[1, 1];
		xlabel="Index",
		ylabel="Relative flux",
	)
	
	comp_names = names.vals[2:4]
	f_wlc_targ = f_target_wlc[:]
	f_wlc_comps = f_comps_wlc ./ mean(f_comps_wlc, dims=1)
	
	for (cName, col) in zip(sort(comp_names), eachcol(f_wlc_comps))
		lines!(ax, col; label=cName)
	end
	
	lines!(ax, f_wlc_targ ./ mean(f_wlc_targ, dims=1);
		linewidth = 5,
		color = :darkgrey,
		label = "WASP-50",
	)
	
	axislegend()
	

	fig #|> as_svg
end

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

# ╔═╡ bd00ebca-4ed2-479b-b1c6-4ba0a7a1043c
md"""
We plot these below for each comparison star division. Move the slider to view the plot for the corresponding comparison star:

comp star $(@bind comp_idx Slider(1:ncomps, show_value=true))

!!! note "Future"
	Ability to interact with sliders completely in the browser coming soon!
"""

# ╔═╡ 49d04cf9-2bec-4350-973e-880e376ab428
begin	
	offs = reshape(range(0, 0.3, length=nbins), 1, :) # Arbitrary offsets for clarity
	f_norm_w = Array{Float64}(undef, ntimes, ncomps, nbins)
	for c_i in 1:ncomps
		f_w = oLCw ./ cLCw[:, c_i, :]
		f_norm_w[:, c_i, :] .= f_w ./ median(f_w, dims=1) .+ offs
	end
	baselines = ones(size(f_norm_w[:, comp_idx, :])) .+ offs # Reference baselines
end;

# ╔═╡ 852d8bdf-96de-4564-83ec-c83d18893970
# Median filtered curves for visualization purposes
f_med, _, f_diff = filt(f_norm_w[:, comp_idx, :], window_width)

# ╔═╡ af8109e7-aac2-480d-a83b-c6dd8c53d380
filt(f_norm_w[:, comp_idx, :], window_width)

# ╔═╡ 97d7166f-bcb2-4ecf-80ed-15d185580cf3
f_norm_w[:, comp_idx, :], window_width

# ╔═╡ 9a6d25a2-6a44-49a7-a9e8-aa3651d67ae0
let
	fig = Figure(resolution=FIG_TALL)
	comp_names = names.vals[begin+1:end]
	ax_left = Axis(fig[1, 1], title = "target / $(comp_names[comp_idx])")
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
	
	# ax.xlabel = "Index"
	# ax.ylabel = "Relative flux + offset"
	
	save(
		"../../ACCESS_WASP-50b/figures/reduced/divided_blcs_LDSS3.pdf",
		fig
	)
	
	fig #|> as_svg
end

# ╔═╡ 5110de9a-3721-4043-b8b7-493daacb4137
target_binned_mags = mapslices(f_to_med_mag, oLCw, dims=1)[use_idxs, :]

# ╔═╡ 53f5a645-93e0-499a-bb36-e4ff0153a63c
comp_binned_mags = mapslices(f_to_med_mag, cLCw, dims=1)[use_idxs, :, :]

# ╔═╡ c03cb527-d16d-47aa-ab63-6970f4ff0b1f
times, airmass = getindex.(
	Ref(LC["temporal"].columns),
	["bjd", "airmass"]
)

# ╔═╡ 354580e4-0aa9-496f-b024-665025a2eeda
lc = let
	med_mag = f_to_med_mag(f_target_wlc |> vec)
	hcat(times, med_mag, zeros(length(times)))[use_idxs, :]
end

# ╔═╡ 898a390b-49f7-45f4-b1a1-b22922d69a29
let
	savepath = "$(tdir)/white-light"
	rm(savepath, force=true, recursive=true)
	mkpath(savepath)
	writedlm("$(savepath)/lc.dat", lc, ",    ")
	writedlm("$(savepath)/comps.dat", comps, ",    ")
end

# ╔═╡ 631c4d02-58cc-4c70-947f-22c8e8b11015
let
	savepath = "$(tdir)/wavelength"
	rm(savepath, force=true, recursive=true)
	for i in 1:nbins
		save_path_w = "$(savepath)/wbin$(i-1)"
		mkpath("$(save_path_w)")
		lc_w = hcat(
			times[use_idxs], target_binned_mags[:, i], zeros(length(times[use_idxs]))
		)
		writedlm("$(save_path_w)/lc.dat", lc_w, ",    ")
		writedlm("$(save_path_w)/comps.dat", comp_binned_mags[:, :, i], ",    ")
	end
end

# ╔═╡ c65c2298-e3a3-4666-be9d-73ee43d94847
fwhm, trace_center, sky_flux = median_eparam.(
	["width", "peak", "sky"],
	Ref(LC["cubes"]),
	Ref(target_name),
	Ref(common_wav_idxs[binned_wav_idxs_range])
)

# ╔═╡ ecc8917c-96a1-488c-9e0b-2e8f6ca6cd60
LC["cubes"]["peak"] |> keys

# ╔═╡ 079c3915-33af-40db-a544-28453732c372
specshifts = load_npz(
	"$(dirname(fpath))/specshifts.npy", allow_pickle=true
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

# ╔═╡ f788835c-8e81-4afe-805e-4caf2d5e5d5b
md"""
## Notebook setup
"""

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
# ╠═34ef4580-bb95-11eb-34c1-25893217f422
# ╟─a8e3b170-3fc2-4d12-b117-07bd37d27710
# ╠═01e27cf6-0a1b-4741-828a-ef8bf7037ae9
# ╠═3ce1d29a-7219-49be-baba-4125b5a557de
# ╠═de51df33-1be5-473b-afa9-c5a2fc40493a
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
# ╠═80480fc2-861b-4fa9-be2d-fcb57f479a93
# ╠═45418bd3-74a3-4758-9fce-adddbeeec076
# ╠═f6ed9c49-bf9c-467d-a602-8d18758a4602
# ╠═2312f321-67ff-4fec-ba73-ff78c82af8dd
# ╠═979895ea-9b13-494c-be08-2496562ccf07
# ╠═ce86c8d1-1900-413f-90be-c9f4eb386a5a
# ╠═bf61788d-ca13-457e-8026-62bc2460d4d7
# ╠═86cedb14-37d1-494b-b022-181887195719
# ╠═80d2f86f-377f-4dac-9b89-89518f001a6a
# ╠═48f31cef-ee18-47f9-a71d-4a5b2c6efe71
# ╠═7d68ad39-3e39-48fa-939a-e56c6659d2b3
# ╠═a6088ea2-904f-4909-b1be-9470e7ec2010
# ╟─bd937d51-17e9-4de3-a5d0-4c436d413940
# ╠═9697e26b-b6d9-413b-869f-47bc2ab99919
# ╠═8c82a78d-7382-40d4-a76b-02e7cd061d67
# ╠═cf38810c-9b0c-4194-bea3-e0aa26e7cf98
# ╠═13385b21-fbd7-484d-a1ac-0687834f92c7
# ╠═cb805821-5d2e-484c-93a5-10a897d2cfe7
# ╠═ded63b4b-61b6-41b6-98d4-d13166bce76a
# ╠═4b2ed9db-0a17-4e52-a04d-3a6a5bf2c054
# ╠═9372c69a-0aad-4e6e-9ea3-e934fa09b758
# ╟─d5c6d058-17c6-4cf0-97b8-d863b1529161
# ╠═20d12d7b-c666-46c3-8f48-5501641e8df3
# ╠═470514e7-0f08-44a3-8519-5d704ea6b8d4
# ╠═f80347e8-dc5a-4b0c-a6c0-db5c12eadcbb
# ╠═fda9022a-2b95-4d06-8a80-80321a1bf441
# ╠═1d91f3b2-293f-467c-ba40-7dc6497152cd
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
# ╠═852d8bdf-96de-4564-83ec-c83d18893970
# ╠═af8109e7-aac2-480d-a83b-c6dd8c53d380
# ╠═97d7166f-bcb2-4ecf-80ed-15d185580cf3
# ╟─bd00ebca-4ed2-479b-b1c6-4ba0a7a1043c
# ╠═9a6d25a2-6a44-49a7-a9e8-aa3651d67ae0
# ╟─af07dc54-eeb5-4fbe-8dd0-289dea97502a
# ╟─88cc640d-b58d-4cde-b793-6c66e74f6b3a
# ╠═301ff07c-8dd5-403a-bae8-a4c38deeb331
# ╠═c03cb527-d16d-47aa-ab63-6970f4ff0b1f
# ╠═c65c2298-e3a3-4666-be9d-73ee43d94847
# ╠═ecc8917c-96a1-488c-9e0b-2e8f6ca6cd60
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
# ╟─3b6b57e2-46ab-46d9-b334-6264daf583f3
# ╠═2d34e125-548e-41bb-a530-ba212c0ca17c
# ╟─f788835c-8e81-4afe-805e-4caf2d5e5d5b
# ╠═c911cecd-0747-4cd1-826f-941f2f58091c
# ╟─b2c61d08-6fcf-4b0c-a21a-c0c5e3205210
