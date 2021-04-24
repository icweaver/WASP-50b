### A Pluto.jl notebook ###
# v0.14.3

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
	import PlutoUI as pl
	using CairoMakie
	using Glob
	using PyCall
	using Statistics
end

# ╔═╡ 2d53ee50-c95b-4e9c-9ddb-e47ed9e3168b
md"""
# LDSS3 reduced data

$(pl.TableOfContents())

In this notebook we will visualize the stellar spectra, white-light, and wavelength binned light curves from the raw flux extracted from LDSS3 using [(cite)](). Select a data cube to begin:

$(@bind data_path pl.Select(sort(glob("data/reduced_LDSS3/*.npy"))))
"""

# ╔═╡ 470c738f-c2c8-4f56-b936-e08c43c161b1
md"""
## Data cube exploration 🔳

This is an `npy` file that contains the following fields:

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
"""

# ╔═╡ 4dddc586-2383-41f8-a888-ef421372d71a
md"""
The contents of each field can also be interactively explored using the tree viewer below:
"""

# ╔═╡ 3653ee36-35a6-4e0a-8d46-4f8389381d45
begin
	py"""
	import numpy as np
	
	def load_npz(fpath, allow_pickle=False):
		return np.load(fpath, allow_pickle=allow_pickle)[()]
	"""
	load_npz(s; allow_pickle=false) = py"load_npz"(s, allow_pickle=allow_pickle)
end

# ╔═╡ 44808b97-df11-4aff-9e97-f97987fe9939
cube = load_npz(data_path, allow_pickle=true)

# ╔═╡ 51898712-b435-476f-8693-6824f77f3af8
data_path

# ╔═╡ 7185f603-3e57-42db-9665-c215094776ad
md"""
Next we extract the common wavelength grid, along with the target and comparison star flux, and compute the resulting spectrum and light curves:
"""

# ╔═╡ 31bdc830-dfc9-445a-ba7e-76be7627561a
wav = cube["spectral"]["wavelength"]

# ╔═╡ 5b645084-1f21-42a2-8184-e27f8b3000c3
target_name = cube["target"]

# ╔═╡ 5d624218-5106-4e3f-a3b1-77ed0d2a02fa
comp_names = convert(Vector{String}, cube["comparisons"])

# ╔═╡ 639f666b-09fc-488b-982b-01a523278cae
raw_counts = cube["cubes"]["raw_counts"];

# ╔═╡ 7c6a3f88-cdc5-4b56-9bbf-2e9a2fa8ae26
target_fluxes = raw_counts[target_name]

# ╔═╡ 7bfc971c-8737-49ad-adec-ac57d176f10e
md"""
We can extract the comparison star flux in a similar way by stacking the ``N \times W`` matrix for each star on top of each other:
"""

# ╔═╡ 968d1f11-12a0-4b8f-abb5-0195652e4e1f
md"""
With the relevant fields extracted, we can now explore the various data products below.
"""

# ╔═╡ e774a20f-2d58-486a-ab71-6bde678b26f8
md"""
## Stellar spectra 🌟
"""

# ╔═╡ ed92fff0-6413-4fc2-939e-ecd13430ce75
md"""
This is the extracted flux from each star, as a function of wavelength.
"""

# ╔═╡ 843fc32b-a0c8-4281-876a-9579a862548b
md"""
### Plot
"""

# ╔═╡ d6c17e99-452e-43e0-b998-10cb81009076
md"""
Move the slider below to view the stellar spectra for WASP-50b and each comparison start at any given time:

!!! future
	Ability to interact with sliders completely in the browser coming soon!
"""

# ╔═╡ e3468c61-782b-4f55-a4a1-9d1883655d11
md"""
## White light curves ⚪
"""

# ╔═╡ 08eafc08-b0fb-4c99-9d4e-7fa5085d386c
md"""
With these spectra now extracted, we can integrate the flux along the spectral direction to build the white-light curves for each star. This will be stored in the `ntimes` ``\times`` 1 ``\times`` `ncomps` array, `f_norm_wlc`.
"""

# ╔═╡ ae39e476-3850-4bbf-aa45-0a16c2425324
md"""
### Integrate spectra
"""

# ╔═╡ e9bbc964-d4e5-4d60-b188-68807e0a030d
md"""
Dividing the target wlc by each comparison star will eliminate some of the common systematics shared between them (e.g., air mass, local refractive atmospheric effects), and make the shape of the transit white-light curve more apparent.
"""

# ╔═╡ 7831034a-01fc-49bb-bd42-ff796713cc50
md"""
### Plot
"""

# ╔═╡ e98dee2e-a369-448e-bfe4-8fea0f318fa8
md"""
## Binned light curves 🌈
"""

# ╔═╡ 83da3243-c61b-4a4b-8d50-6d32c606d34c
md"""
Similarly to above, we next integrate the target and comparison star flux over given wavelength ranges to build the binned light curves below.
"""

# ╔═╡ e992642a-342b-467b-ad75-594b73c88ac7
md"""
### Define wavelength bins
"""

# ╔═╡ 17039dc4-d554-4810-a147-854cd95c81d3
md"""
We define the lower and upper bound for each of the `nbins` wavelength bins that the timeseries data will be integrated over, where `nbins` is the total number of wavelength bins used:
"""

# ╔═╡ 84893ece-b867-4e13-a45f-081d4a9a0934
const λ_start, Δλ, nbins = 5_800, 200, 19

# ╔═╡ 9aa02a13-ec25-40cc-8a54-89a192d228e5
wbin_idx_start = findfirst(==(λ_start), wav)

# ╔═╡ 1b550fd4-95a8-4130-972f-8dc963436cdc
compute_wbins_idxs(c, λ_start_idx, Δλ) = (λ_start_idx:λ_start_idx+Δλ) .+ c

# ╔═╡ fadb7a1b-d631-493f-97e0-4271599cd480
wbins_idxs = compute_wbins_idxs.(0:Δλ:Δλ*nbins-1, wbin_idx_start, Δλ)

# ╔═╡ a0b2cff5-67cd-4ee0-b538-d39fe6c6b67c
wav[wbins_idxs[begin][begin]], wav[wbins_idxs[end][end]]

# ╔═╡ f7feb44e-a363-4f8d-bf62-d3541533e4da
md"""
### Compute binned divided LCs
"""

# ╔═╡ eb4f7a92-a9e7-4bf5-8b1c-bca928fcedde
md"""
We will next compute `oLCw` and `cLCw`, where `oLCw` is an `ntimes` ``\times`` `nbins` matrix that holds the binned target flux, where `ntimes` is the number of timeseries points ``N``. Similarly `cLCw` is an `ntimes` ``\times`` `nbins` ``\times`` `ncomps` matrix that holds the comparison star flux, where `ncomps` is the number of comparison stars:
"""

# ╔═╡ 6953f995-86e1-42ee-9591-8c4e372b06fb
const ntimes = size(target_fluxes)[1]

# ╔═╡ 081c0175-ea02-4a41-9dfe-60b871efc27d
md"""
Time index: $(@bind t_i pl.Slider(1:10:ntimes, show_value=true))
"""

# ╔═╡ 70a2863d-6f55-4cec-8c69-162157e2c199
const ncomps = length(comp_names)

# ╔═╡ bb2c6085-5cb4-4efc-a322-27fda09fb904
comp_fluxes = cat([raw_counts[comp_names[c_i]] for c_i in 1:ncomps]..., dims=3)

# ╔═╡ 2fd7bc68-a6ec-4ccb-ad73-07b031ffef5a
let
	fig = Figure(resolution = (1200, 700))

	lines(fig[1, 1], wav, target_fluxes[t_i, :], label="WASP-50")
	lines(fig[2, 1], wav, comp_fluxes[t_i, :, 1], label="comp 1")
	lines(fig[1, 2], wav, comp_fluxes[t_i, :, 2], label="comp 2")
	lines(fig[2, 2], wav, comp_fluxes[t_i, :, 3], label="comp 3")
	
	axs = reshape(copy(fig.content), 2, 2)
	axislegend.(axs)
	
	scene = current_axis()
	linkaxes!(axs...)
	hidexdecorations!.(axs[1, :], grid=false)
	hideydecorations!.(axs[:, 2], grid=false)
	ylims!(scene, 0, 8e4)
	
	fig[1:2, 0] = Label(fig, "Counts", rotation=π/2)
	fig[end+1, 2:3] = Label(fig, "Index")
    #textsize = 30, font = noto_sans_bold, color = (:black, 0.25))
	
	fig
end

# ╔═╡ 756c0e88-393e-404f-b36c-9f861a8172e2
begin
	f_wlc = sum(target_fluxes, dims=2) ./ sum(comp_fluxes, dims=2)
	f_norm_wlc = f_wlc ./ median(f_wlc, dims=1) 
end;

# ╔═╡ 76ba54ff-a396-48e2-8837-b839cb62caef
let
	fig = Figure(resolution = (1200, 1_000))
	
	for c_i in 1:ncomps
		scatter(
			fig[c_i, 1], vcat(f_norm_wlc[:, :, c_i]...),
			strokewidth=0,
			label="target / comp $(c_i)"
		)
	end
	
	axs = copy(fig.content)
	axislegend.(axs, valign=:top)
	
	scene = current_axis()
	linkaxes!(axs...)
	hidexdecorations!.(axs[1:2], grid = false)
	axs[end].xlabel = "Index"
	axs[2].ylabel = "Relative flux"
	ylims!(scene, 0.97, 1.02)
	
	fig
end

# ╔═╡ 2b20e471-16e5-4b54-abc5-4308af4e60b6
"""
	sum_flux(A::AbstractMatrix, idx_range::AbstractRange, [dims=2])

For columnar data `A`, returns the sum of each row (for default `dims=2`) along the range of columns specified by `idx_range`.
"""
sum_flux(A, idx_range, dims=2) = sum(view(A, :, idx_range), dims=dims)[:, 1]

# ╔═╡ a6805e63-bbf4-48bc-8e15-be24da9b0348
begin
	oLCw = Matrix{Float64}(undef, ntimes, nbins)
	cLCw = Array{Float64, 3}(undef, ntimes, nbins, ncomps)
	
	for w_i in 1:nbins
		oLCw[:, w_i] .= sum_flux(target_fluxes, wbins_idxs[w_i])
		for c_i in 1:ncomps
			cLCw[:, w_i, c_i] .= sum_flux(comp_fluxes[:, :, c_i], wbins_idxs[w_i])
		end
	end
end

# ╔═╡ f44393ab-2e49-4817-9870-890c9cd556ce
md"""
With `oLCw` and `cLCw` now computed, we next compute `f_norm_w`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\times`` `nbins` ``\times`` `ncomps`, where, for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve. We plot these below for each comparison star division:
"""

# ╔═╡ e59f0d21-7573-43c6-9eb9-8863f174f77e
md"""
### Plot
"""

# ╔═╡ 70380521-e0f3-45a4-b818-c50adb635c69
md"""
Move the slider to view the plot for the corresponding comparison star:

comp star $(@bind comp_idx pl.Slider(1:ncomps, show_value=true))

!!! future
	Ability to interact with sliders completely in the browser coming soon!
"""

# ╔═╡ 90fec9cf-37b0-4bcf-a6df-85a4cdfc511b
begin	
	f_norm_w = Array{Float64}(undef, ntimes, nbins, ncomps)
	for c_i in 1:ncomps
		f_w = oLCw ./ cLCw[:, :, c_i]
		f_norm_w[:, :, c_i] .= f_w ./ median(f_w, dims=1)
	end
	
	offs = reshape(range(0, 0.6, length=nbins), 1, :) # Arbitrary offsets for clarity
	f_norm_w .+= offs
	baselines = ones(size(f_norm_w[:, :, comp_idx])) .+ offs # Reference baselines
end;

# ╔═╡ fbc57d8b-3b1b-44d1-bd7d-0e9749026d4c
begin
	fig = Figure(resolution = (600, 900))
	ax = fig[1, 1] = Axis(fig, title = "target / comp $comp_idx")
	
	cmap = reverse(to_colormap(:Spectral_4, nbins))
	
	for (c, f, b) in zip(cmap, eachcol(f_norm_w[:, :, comp_idx]), eachcol(baselines))
		scatter!(f, color=c, strokewidth=0, markersize=5)
		lines!(b, color=c)
	end
		
	scene = current_axis()
	ylims!(scene, 0.9, 1.7)
	ax.xlabel = "Index"
	ax.ylabel = "Relative flux + offset"
	
	fig
end

# ╔═╡ Cell order:
# ╟─2d53ee50-c95b-4e9c-9ddb-e47ed9e3168b
# ╟─470c738f-c2c8-4f56-b936-e08c43c161b1
# ╟─4dddc586-2383-41f8-a888-ef421372d71a
# ╠═44808b97-df11-4aff-9e97-f97987fe9939
# ╠═3653ee36-35a6-4e0a-8d46-4f8389381d45
# ╠═51898712-b435-476f-8693-6824f77f3af8
# ╟─7185f603-3e57-42db-9665-c215094776ad
# ╠═31bdc830-dfc9-445a-ba7e-76be7627561a
# ╠═5b645084-1f21-42a2-8184-e27f8b3000c3
# ╠═5d624218-5106-4e3f-a3b1-77ed0d2a02fa
# ╠═639f666b-09fc-488b-982b-01a523278cae
# ╠═7c6a3f88-cdc5-4b56-9bbf-2e9a2fa8ae26
# ╟─7bfc971c-8737-49ad-adec-ac57d176f10e
# ╠═bb2c6085-5cb4-4efc-a322-27fda09fb904
# ╟─968d1f11-12a0-4b8f-abb5-0195652e4e1f
# ╟─e774a20f-2d58-486a-ab71-6bde678b26f8
# ╟─ed92fff0-6413-4fc2-939e-ecd13430ce75
# ╟─843fc32b-a0c8-4281-876a-9579a862548b
# ╟─d6c17e99-452e-43e0-b998-10cb81009076
# ╟─081c0175-ea02-4a41-9dfe-60b871efc27d
# ╠═2fd7bc68-a6ec-4ccb-ad73-07b031ffef5a
# ╟─e3468c61-782b-4f55-a4a1-9d1883655d11
# ╟─08eafc08-b0fb-4c99-9d4e-7fa5085d386c
# ╟─ae39e476-3850-4bbf-aa45-0a16c2425324
# ╠═756c0e88-393e-404f-b36c-9f861a8172e2
# ╟─e9bbc964-d4e5-4d60-b188-68807e0a030d
# ╟─7831034a-01fc-49bb-bd42-ff796713cc50
# ╠═76ba54ff-a396-48e2-8837-b839cb62caef
# ╟─e98dee2e-a369-448e-bfe4-8fea0f318fa8
# ╟─83da3243-c61b-4a4b-8d50-6d32c606d34c
# ╟─e992642a-342b-467b-ad75-594b73c88ac7
# ╟─17039dc4-d554-4810-a147-854cd95c81d3
# ╠═84893ece-b867-4e13-a45f-081d4a9a0934
# ╠═a0b2cff5-67cd-4ee0-b538-d39fe6c6b67c
# ╠═9aa02a13-ec25-40cc-8a54-89a192d228e5
# ╠═1b550fd4-95a8-4130-972f-8dc963436cdc
# ╠═fadb7a1b-d631-493f-97e0-4271599cd480
# ╟─f7feb44e-a363-4f8d-bf62-d3541533e4da
# ╟─eb4f7a92-a9e7-4bf5-8b1c-bca928fcedde
# ╠═6953f995-86e1-42ee-9591-8c4e372b06fb
# ╠═70a2863d-6f55-4cec-8c69-162157e2c199
# ╠═2b20e471-16e5-4b54-abc5-4308af4e60b6
# ╠═a6805e63-bbf4-48bc-8e15-be24da9b0348
# ╟─f44393ab-2e49-4817-9870-890c9cd556ce
# ╠═90fec9cf-37b0-4bcf-a6df-85a4cdfc511b
# ╟─e59f0d21-7573-43c6-9eb9-8863f174f77e
# ╟─70380521-e0f3-45a4-b818-c50adb635c69
# ╠═fbc57d8b-3b1b-44d1-bd7d-0e9749026d4c
# ╠═b1b0690a-a1eb-11eb-1590-396d92c80c23
