### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 691eddff-f2eb-41a8-ab05-63afb46d15f2
begin
	import PlutoUI as pl
	using CairoMakie
	using Glob
	using PyCall
	using Statistics, KernelDensity
end

# ╔═╡ 506eeeb2-e56d-436b-91b8-605e52201563
md"""
# Detrended white light curves

Gonna do some stuff

$(pl.TableOfContents())
"""

# ╔═╡ a8cf11e2-796e-45ff-bdc9-e273b927700e
md"""
## Transit curves ⚪
"""

# ╔═╡ 39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
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
end

# ╔═╡ 2191791b-df62-4f1b-88bf-060cc47896b2
cube = load_npz(
	"data/detrended_wlcs/out_l/WASP50/w50_161211/white-light/BMA_WLC.npy",
	allow_pickle = true
)

# ╔═╡ 4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
begin
	scatter(cube["t"] .- 2.45e6, cube["LC_det"])
	lines!(cube["t_interp"] .- 2.45e6, cube["LC_det_model_interp"])
	current_figure()
end

# ╔═╡ 68ec4343-5f6c-4dfd-90b5-6393b4c819b9
md"""
## Corner plots 📐
"""

# ╔═╡ 931ce3d5-c4ed-496c-883b-d7ee33e957cc
d = filter!(
 	p -> p.first ∈ ["p", "t0", "P", "rho", "inc", "b", "aRs", "q1"],
	load_pickle(
		"data/detrended_wlcs/out_l/WASP50/w50_131219/white-light/BMA_posteriors.pkl"
	)
)

# ╔═╡ c4524acb-4656-47a4-850f-f8ff1408b435
pair = cat(d["rho"], d["inc"], dims=2)

# ╔═╡ f65babf8-7d6a-4528-b33b-1d71d2047cc6
begin
	fig = Figure()
	
	hist(fig[1, 1], pair[:, 1])
	contourf(fig[2, 1], kde(pair, npoints=(16, 16)), levels=10)
	hist(fig[2, 2], pair[:, 2])
	
	fig
end

# ╔═╡ baeadfce-535a-46c3-8cb9-79cf6bde8555
md"""
## Packages 📦
"""

# ╔═╡ Cell order:
# ╟─506eeeb2-e56d-436b-91b8-605e52201563
# ╠═a8cf11e2-796e-45ff-bdc9-e273b927700e
# ╠═39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
# ╠═2191791b-df62-4f1b-88bf-060cc47896b2
# ╠═4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
# ╟─68ec4343-5f6c-4dfd-90b5-6393b4c819b9
# ╠═931ce3d5-c4ed-496c-883b-d7ee33e957cc
# ╠═c4524acb-4656-47a4-850f-f8ff1408b435
# ╠═f65babf8-7d6a-4528-b33b-1d71d2047cc6
# ╟─baeadfce-535a-46c3-8cb9-79cf6bde8555
# ╠═691eddff-f2eb-41a8-ab05-63afb46d15f2
