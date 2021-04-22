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
	import PlutoUI
	#using Gadfly, ColorSchemes, Compose
	#import Cairo, Fontconfig
	using Statistics
	using PyCall
end

# ╔═╡ 4a0b5c46-aebd-4f72-8873-df7094c45c78
using ColorSchemes

# ╔═╡ 4962f384-5882-4932-a487-6a5a3b7b0ef6
begin
	using CairoMakie
end

# ╔═╡ 2d53ee50-c95b-4e9c-9ddb-e47ed9e3168b
md"""
# Extracted WLCs
"""

# ╔═╡ 1293b7f7-3054-4504-965f-2bd8ebe05e9e
TableOfContents()

# ╔═╡ 470c738f-c2c8-4f56-b936-e08c43c161b1
md"""
## LDSS3 data cube exploration

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
λ_start, Δλ, nbins = 5_800, 199, 20

# ╔═╡ 1b550fd4-95a8-4130-972f-8dc963436cdc
compute_wbins_idxs(c, λ_start_idx, Δλ) = (λ_start_idx:λ_start_idx+Δλ) .+ c

# ╔═╡ f7feb44e-a363-4f8d-bf62-d3541533e4da
md"""
### Compute binned, divided LCs
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
With `oLCw` and `cLCw` now computed, we next compute `f_norms`, the binned target flux divided by each comparison star binned flux, normalized by the median of the original ratio. This has dimensions `ntimes` ``\times`` `nbins` ``\times`` `ncomps`, where, for a given comparison star, each column from the resulting matrix corresponds to a final binned light curve:
"""

# ╔═╡ f3644805-5951-4e5d-8d20-1596eebcce65
# begin
# 	wbin_labels = [
# 		"$(wav[w_idxs[begin]]) - $(wav[w_idxs[end]]) Å" for w_idxs in wbins_idxs
# 	]
	
# 	offs = reshape(range(0, 0.4, length=nbins), 1, :)
	
# 	p_blcs = plot(
# 		# Binned LCs
# 		layer(
# 			f_norms[:, :, comp_idx] .+ offs,
# 			x = Row.index,
# 			y = Col.value,
# 			color = Col.index,
# 		),
# 		# Baseline
# 		layer(
# 			ones(size(f_norms[:, :, comp_idx])).+ offs,
# 			x = Row.index,
# 			y = Col.value,
# 			color = Col.index,
# 			Geom.line,
# 		),
# 		# Label wavelength bins
# 		Guide.annotation(
# 			compose(
# 				context(),
# 				text(repeat([ntimes], nbins), offs .+ 1.01, wbin_labels, [hright])
# 			)	
# 		),
# 		# Global settings
# 		Theme(highlight_width=0pt, key_position=:none),
# 		Scale.ContinuousColorScale(
# 			p -> get(reverse(ColorSchemes.Spectral_4), p)
# 		),
# 		Guide.xlabel("Index"),
# 		Guide.ylabel("Normalized flux + offset"),
# 		Guide.title("Reduced target flux / comp $comp_idx flux"),
# 	) #|> HTMLDocument
	
# 	#plot(y=ones(size(f_norms[:, :, 1])).+ offs)
# 	#draw(PNG("/home/mango/Desktop/myplot.png", 8inch, 11inch), p_blcs)
# end

# ╔═╡ 3653ee36-35a6-4e0a-8d46-4f8389381d45
begin
	py"""
	import numpy as np
	import pickle
	
	def load_pickle(fpath):
		with open(fpath, "rb") as f:
			data = pickle.load(f)
		return data
	
	def load_npz(fpath, allow_pickle=False):
		return np.load(fpath, allow_pickle=allow_pickle)[()]
	"""
	load_pickle(s) = py"load_pickle"(s)
	load_npz(s; allow_pickle=false) = py"load_npz"(s, allow_pickle=allow_pickle)
end

# ╔═╡ 44808b97-df11-4aff-9e97-f97987fe9939
cube = load_npz(
	"/home/mango/data/WASP50/LDSS3/spectralCube_WASP-50b_LDSS3_e150927_4stars_316spectra_05px_shifted.npy",
	allow_pickle = true,
)

# ╔═╡ 5b645084-1f21-42a2-8184-e27f8b3000c3
target_name = cube["target"]

# ╔═╡ 1f3f096f-0aa3-44bc-b52b-568da4e75cf3
comp_names = convert(Vector{String}, cube["comparisons"])

# ╔═╡ 70a2863d-6f55-4cec-8c69-162157e2c199
ncomps = length(comp_names)

# ╔═╡ 6fb2ce35-23ab-4a14-85b1-b8c4fc53d15d
md"""
We plot these below for a given comparison star, which can be selected from slider: $(@bind comp_idx Slider(1:ncomps, show_value=true) )
"""

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

# ╔═╡ 6953f995-86e1-42ee-9591-8c4e372b06fb
ntimes = size(target_fluxes)[1]

# ╔═╡ bb2c6085-5cb4-4efc-a322-27fda09fb904
comp_fluxes = [raw_counts[comp_names[c_i]] for c_i in 1:ncomps]

# ╔═╡ a6805e63-bbf4-48bc-8e15-be24da9b0348
begin
	oLCw = Matrix{Float64}(undef, ntimes, nbins)
	cLCw = Array{Float64, 3}(undef, ntimes, nbins, ncomps)
	
	for w_i in 1:nbins
		oLCw[:, w_i] .= sum_flux(target_fluxes, wbins_idxs[w_i])
		for c_i in 1:ncomps
			cLCw[:, w_i, c_i] .= sum_flux(comp_fluxes[c_i], wbins_idxs[w_i])
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

# ╔═╡ ef7eade9-5c9b-4d94-8b18-5ad208729e7f
begin
	wbin_labels = [
		"$(wav[w_idxs[begin]]) - $(wav[w_idxs[end]]) Å" for w_idxs in wbins_idxs
	]
	
	offs = reshape(range(0, 0.4, length=nbins), 1, :)
	lines(f_norms[:, :, comp_idx] .+ offs)
		
	current_figure()
end

# ╔═╡ 7ee58880-d19e-4ecc-8ad2-74b057c8ebdc
struct HTMLDocument
	embedded
end

# ╔═╡ 6199f3e8-a082-4b36-bcae-44a3a4e998ee
function Base.show(io::IO, mime::MIME"text/html", doc::HTMLDocument)
	println(io, "<html>")
	show(io, mime, doc.embedded)
	println(io, "</html>")
end

# ╔═╡ efccdc2d-0ed3-47ee-9a3d-d7bfcfde3fa8
# begin
# 	using JSServe
# 	Page()
# end

# ╔═╡ cc6547eb-548a-4f4a-9081-32a19d52acb2
X = [
	6 7 8 9
	1 2 4 8
	9 8 9 7
	6 6 6 6
	1 1 1 1
	0 3 4 1
]

# ╔═╡ ad379c1d-4efe-48b4-89ef-ed03cbbdbce6
function CairoMakie.lines(A::Matrix)
	nrows, ncols = size(A)
	colors = range(ColorSchemes.Spectral_4[1], ColorSchemes.Spectral_4[4], length=ncols)
	a, l = lines([0])
	for (line, color) in zip(eachcol(A), colors)
		lines!(l, line, color=color)
	end
	ylims!(l, 0.8, 1.6)
	#axislegend()
	#current_figure()
end

# ╔═╡ b5e8d9d2-14d6-4d92-b051-837f24daf335
begin
	lines([1,2,3])
	ylims!(current_figure(), (0, 1))
end

# ╔═╡ Cell order:
# ╟─2d53ee50-c95b-4e9c-9ddb-e47ed9e3168b
# ╟─1293b7f7-3054-4504-965f-2bd8ebe05e9e
# ╟─470c738f-c2c8-4f56-b936-e08c43c161b1
# ╠═44808b97-df11-4aff-9e97-f97987fe9939
# ╠═5b645084-1f21-42a2-8184-e27f8b3000c3
# ╠═1f3f096f-0aa3-44bc-b52b-568da4e75cf3
# ╠═31bdc830-dfc9-445a-ba7e-76be7627561a
# ╠═639f666b-09fc-488b-982b-01a523278cae
# ╠═7c6a3f88-cdc5-4b56-9bbf-2e9a2fa8ae26
# ╠═bb2c6085-5cb4-4efc-a322-27fda09fb904
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
# ╟─6fb2ce35-23ab-4a14-85b1-b8c4fc53d15d
# ╠═ef7eade9-5c9b-4d94-8b18-5ad208729e7f
# ╠═f3644805-5951-4e5d-8d20-1596eebcce65
# ╠═3653ee36-35a6-4e0a-8d46-4f8389381d45
# ╠═7ee58880-d19e-4ecc-8ad2-74b057c8ebdc
# ╠═6199f3e8-a082-4b36-bcae-44a3a4e998ee
# ╠═b1b0690a-a1eb-11eb-1590-396d92c80c23
# ╠═efccdc2d-0ed3-47ee-9a3d-d7bfcfde3fa8
# ╠═cc6547eb-548a-4f4a-9081-32a19d52acb2
# ╠═ad379c1d-4efe-48b4-89ef-ed03cbbdbce6
# ╠═b5e8d9d2-14d6-4d92-b051-837f24daf335
# ╠═4a0b5c46-aebd-4f72-8873-df7094c45c78
# ╠═4962f384-5882-4932-a487-6a5a3b7b0ef6
