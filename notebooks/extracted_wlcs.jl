### A Pluto.jl notebook ###
# v0.14.2

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
	using PlutoUI
	using Statistics#, StatsPlots
	using PyCall
end

# ╔═╡ 36e5bb5f-75a1-47c3-baa3-66bd77eb9515
using Gadfly

# ╔═╡ 2d53ee50-c95b-4e9c-9ddb-e47ed9e3168b
md"""
# LDSS3 reduced data
"""

# ╔═╡ 1293b7f7-3054-4504-965f-2bd8ebe05e9e
md"""
$(TableOfContents())

In this notebook we will visualize the white-light and wavelength binned, divided light curves from the raw flux extracted from LDSS3 using [(cite)]().
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
        * `fractionofapixel`: **Same length array, but not sure what each value here represents?**
* comparisons
    * Similar to `stellar`, but a list of the comparison star square names only
* target
    * Similar to `stellar`, but a list of the target star square name(s) only
* cubes
    * Dict with keys: `raw_counts`, `peak`, `ok`, `sky`, `width`, and `centroid`. Each key contains an ``N \times W`` matrix of measurements, 1 matrix for each square. **I was wondering if the definition for each key is stored somewhere online?**

!!! note
	We will use `wavelength` and `raw_counts` for our analysis here, where `raw_counts` is the sky subtracted, cross-dispersion integrated flux from each star.))
"""

# ╔═╡ 7185f603-3e57-42db-9665-c215094776ad
md"""
Let's extract the target and comparison star flux and compute the resulting spectrum and light curves:
"""

# ╔═╡ 7bfc971c-8737-49ad-adec-ac57d176f10e
md"""
We can extract the comparison star flux in a similar way by stacking the ``N \times W`` matrix for each star on top of each other:
"""

# ╔═╡ 968d1f11-12a0-4b8f-abb5-0195652e4e1f
md"""
With the releavnt fields extracted, we can now explore the various data products below:
"""

# ╔═╡ e774a20f-2d58-486a-ab71-6bde678b26f8
md"""
## Stellar spectra 🌟
"""

# ╔═╡ d6c17e99-452e-43e0-b998-10cb81009076
md"""
Move the slider below to view the stellar spectra for WASP-50b and each comparison start at any given time:
"""

# ╔═╡ e3468c61-782b-4f55-a4a1-9d1883655d11
md"""
## White light curves ⚪
"""

# ╔═╡ 1632c9d9-d61f-4bf6-bbfa-4d5f7deee5f3
let
	p1 = plot(rand(100))
	p2 = plot(rand(100))
	
	testplot = plot(p1, p2, link=:both, ylabel="Y", layout=(2,1))
end

# ╔═╡ e98dee2e-a369-448e-bfe4-8fea0f318fa8
md"""
## Binned light curves 🌈
"""

# ╔═╡ e992642a-342b-467b-ad75-594b73c88ac7
md"""
### Define wavelength bins
"""

# ╔═╡ 17039dc4-d554-4810-a147-854cd95c81d3
md"""
We define the lower and upper bound for each of the `nbins` wavelength bins that the timeseries data will be integrated over:
"""

# ╔═╡ 84893ece-b867-4e13-a45f-081d4a9a0934
const λ_start, Δλ, nbins = 5_800, 200, 19

# ╔═╡ 1b550fd4-95a8-4130-972f-8dc963436cdc
compute_wbins_idxs(c, λ_start_idx, Δλ) = (λ_start_idx:λ_start_idx+Δλ) .+ c

# ╔═╡ f7feb44e-a363-4f8d-bf62-d3541533e4da
md"""
### Compute binned divided LCs
"""

# ╔═╡ eb4f7a92-a9e7-4bf5-8b1c-bca928fcedde
md"""
We will next compute `oLCw` and `cLCw`, where `oLCw` is an `ntimes` ``\times`` `nbins` matrix that holds the binned target flux, where `ntimes` is the number of timeseries points ``N``. Similarly `cLCw` is an `ntimes` ``\times`` `nbins` ``\times`` `ncomps` matrix that holds the comparison star flux, where `ncomps` is the number of comparison stars:
"""

# ╔═╡ c7792d3e-52f9-4762-8e51-00d29d196d93
view_cols(A, col_range) = @view A[:, col_range]

# ╔═╡ 2b20e471-16e5-4b54-abc5-4308af4e60b6
sum_flux(A, idx_range) = sum(view_cols(A, idx_range), dims=2)[:, 1]

# ╔═╡ f44393ab-2e49-4817-9870-890c9cd556ce
md"""
With `oLCw` and `cLCw` now computed, we next compute `f_norms`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\times`` `nbins` ``\times`` `ncomps`, where, for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve. We plot these below for each comparison star division:
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
cube = load_npz("data/reduced_LDSS3/ut150927_flat.npy", allow_pickle = true)

# ╔═╡ 5b645084-1f21-42a2-8184-e27f8b3000c3
target_name = cube["target"]

# ╔═╡ 1f3f096f-0aa3-44bc-b52b-568da4e75cf3
comp_names = convert(Vector{String}, cube["comparisons"])

# ╔═╡ 70a2863d-6f55-4cec-8c69-162157e2c199
const ncomps = length(comp_names)

# ╔═╡ 31bdc830-dfc9-445a-ba7e-76be7627561a
wav = cube["spectral"]["wavelength"]

# ╔═╡ 9aa02a13-ec25-40cc-8a54-89a192d228e5
wbin_idx_start = findfirst(==(λ_start), wav)

# ╔═╡ fadb7a1b-d631-493f-97e0-4271599cd480
wbins_idxs = compute_wbins_idxs.(0:Δλ:Δλ*nbins-1, wbin_idx_start, Δλ)

# ╔═╡ a0b2cff5-67cd-4ee0-b538-d39fe6c6b67c
wav[wbins_idxs[begin][begin]], wav[wbins_idxs[end][end]]

# ╔═╡ 639f666b-09fc-488b-982b-01a523278cae
raw_counts = cube["cubes"]["raw_counts"]

# ╔═╡ 7c6a3f88-cdc5-4b56-9bbf-2e9a2fa8ae26
target_fluxes = raw_counts[target_name]

# ╔═╡ 756c0e88-393e-404f-b36c-9f861a8172e2
wlc_targ = sum(target_fluxes, dims=2)

# ╔═╡ 6953f995-86e1-42ee-9591-8c4e372b06fb
const ntimes = size(target_fluxes)[1]

# ╔═╡ 081c0175-ea02-4a41-9dfe-60b871efc27d
md"""
Time index: $(@bind t_i Slider(1:10:ntimes, show_value=true))
"""

# ╔═╡ bb2c6085-5cb4-4efc-a322-27fda09fb904
comp_fluxes = cat([raw_counts[comp_names[c_i]] for c_i in 1:ncomps]..., dims=3)

# ╔═╡ 2fd7bc68-a6ec-4ccb-ad73-07b031ffef5a
begin
	fig_targ = plot(x=wav, y=target_fluxes[t_i, :], Geom.line)
	fig_comp_1 = plot(x=wav, y=comp_fluxes[t_i, :, 1], Geom.line)
	fig_comp_2 = plot(x=wav, y=comp_fluxes[t_i, :, 2], Geom.line)
	fig_comp_3 = plot(x=wav, y=comp_fluxes[t_i, :, 3], Geom.line)
	gridstack([fig_targ fig_comp_1; fig_comp_2 fig_comp_3])
end

# ╔═╡ 7eec7f0e-0faa-42bd-ae70-23d331604c7d
begin
	p_targ = plot(wav, target_fluxes[t_i, :], label="WASP-50b")
	p_comp_1 = plot(wav, comp_fluxes[t_i, :, 1], label="comp 1")
	p_comp_2 = plot(wav, comp_fluxes[t_i, :, 2], label="comp 2")
	p_comp_3 = plot(wav, comp_fluxes[t_i, :, 3], label="comp 3")
	plot(p_targ, p_comp_1, p_comp_2, p_comp_3, layout=(2, 2), link=:all)
end

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

# ╔═╡ 90fec9cf-37b0-4bcf-a6df-85a4cdfc511b
begin	
	f_norms = Array{Float64}(undef, ntimes, nbins, ncomps)
	for c_i in 1:ncomps
		f = oLCw ./ cLCw[:, :, c_i]
		f_norms[:, :, c_i] .= f ./ median(f, dims=1)
	end
end;

# ╔═╡ 26f86aec-1a39-42f0-9d85-70efe0a490a3
function plot_bins(comp_idx, offs, pal, anns)	
	p = plot(
		title = "target / comp $comp_idx",
		#xlabel = "Index",
		#ylabel = "Relative flux + offset",
	)
	scatter!(p, f_norms[:, :, comp_idx] .+ offs, msw=0, label=false, palette=pal)
	plot!(p, ones(size(f_norms[:, :, comp_idx])) .+ offs, label=false)
	annotate!(p, anns)
	return p
end

# ╔═╡ 981e588f-bb45-4c35-93a3-e4d13ee0181e
begin
	# Define offset between light curves for clarity
	offs = reshape(range(0, 0.6, length=nbins), 1, :)
	
	# Wavelength bin labels
	wbin_labels = [
	"$(wav[w_idxs[begin]]) - $(wav[w_idxs[end]]) Å" for w_idxs in wbins_idxs
	]
	anns = Vector{Tuple{Float64, Float64, Plots.PlotText}}()
	for (wbin_label, off) in zip(wbin_labels, offs)
		push!(
			anns,
			(ntimes, off + 1.01, text("$wbin_label", 10, :right, :bottom))
		)
	end
	
	# Color palette
	pal = palette(:Spectral_4, nbins, rev=:true)
	
	# Plot each binned LC, one comp star per column
	p1 = plot_bins(1, offs, pal, anns)
	p2 = plot_bins(2, offs, pal, anns)
	p3 = plot_bins(3, offs, pal, anns)
	plot(p1, p2, p3, layout=(1, 3), link=:all, size=(1_000, 1_000))
end

# ╔═╡ Cell order:
# ╟─2d53ee50-c95b-4e9c-9ddb-e47ed9e3168b
# ╟─1293b7f7-3054-4504-965f-2bd8ebe05e9e
# ╟─470c738f-c2c8-4f56-b936-e08c43c161b1
# ╠═44808b97-df11-4aff-9e97-f97987fe9939
# ╟─7185f603-3e57-42db-9665-c215094776ad
# ╠═5b645084-1f21-42a2-8184-e27f8b3000c3
# ╠═1f3f096f-0aa3-44bc-b52b-568da4e75cf3
# ╠═31bdc830-dfc9-445a-ba7e-76be7627561a
# ╠═639f666b-09fc-488b-982b-01a523278cae
# ╠═7c6a3f88-cdc5-4b56-9bbf-2e9a2fa8ae26
# ╟─7bfc971c-8737-49ad-adec-ac57d176f10e
# ╠═bb2c6085-5cb4-4efc-a322-27fda09fb904
# ╟─968d1f11-12a0-4b8f-abb5-0195652e4e1f
# ╟─e774a20f-2d58-486a-ab71-6bde678b26f8
# ╟─d6c17e99-452e-43e0-b998-10cb81009076
# ╟─081c0175-ea02-4a41-9dfe-60b871efc27d
# ╠═2fd7bc68-a6ec-4ccb-ad73-07b031ffef5a
# ╠═7eec7f0e-0faa-42bd-ae70-23d331604c7d
# ╟─e3468c61-782b-4f55-a4a1-9d1883655d11
# ╠═756c0e88-393e-404f-b36c-9f861a8172e2
# ╠═1632c9d9-d61f-4bf6-bbfa-4d5f7deee5f3
# ╟─e98dee2e-a369-448e-bfe4-8fea0f318fa8
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
# ╟─c7792d3e-52f9-4762-8e51-00d29d196d93
# ╠═2b20e471-16e5-4b54-abc5-4308af4e60b6
# ╠═a6805e63-bbf4-48bc-8e15-be24da9b0348
# ╟─f44393ab-2e49-4817-9870-890c9cd556ce
# ╠═90fec9cf-37b0-4bcf-a6df-85a4cdfc511b
# ╠═981e588f-bb45-4c35-93a3-e4d13ee0181e
# ╠═26f86aec-1a39-42f0-9d85-70efe0a490a3
# ╠═3653ee36-35a6-4e0a-8d46-4f8389381d45
# ╠═b1b0690a-a1eb-11eb-1590-396d92c80c23
# ╠═36e5bb5f-75a1-47c3-baa3-66bd77eb9515
