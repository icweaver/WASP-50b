### A Pluto.jl notebook ###
# v0.14.3

using Markdown
using InteractiveUtils

# ╔═╡ 1022b4f3-92cb-4b23-a51a-fdd0ef2ff730
begin
	import PlutoUI as pl
	using CairoMakie
	using Glob
	using PyCall
	using Statistics
end

# ╔═╡ 39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
begin
	py"""
	import numpy as np
	
	def load_npz(fpath, allow_pickle=False):
		return np.load(fpath, allow_pickle=allow_pickle)[()]
	"""
	load_npz(s; allow_pickle=false) = py"load_npz"(s, allow_pickle=allow_pickle)
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

# ╔═╡ Cell order:
# ╠═39dbca86-a4b9-11eb-1c64-9ddf1a9990ab
# ╠═2191791b-df62-4f1b-88bf-060cc47896b2
# ╠═4be0d7b7-2ea5-4c4d-92b9-1f8109014e12
# ╠═1022b4f3-92cb-4b23-a51a-fdd0ef2ff730
