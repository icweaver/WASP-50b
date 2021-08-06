### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 239a91a6-f68a-11eb-14fd-0ba8d08b08f9
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using CCDReduction: fitscollection
	using Colors
	using DataFrames
	using DataFramesMeta
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using KernelDensity
	using Latexify
	using LaTeXStrings
	using Measurements
	using Measurements: value, uncertainty
	using NaturalSort
	using OrderedCollections
	using Printf
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
	using Unitful
	
	# Python setup
	ENV["PYTHON"] = "/home/mango/miniconda3/envs/WASP50/bin/python"
	Pkg.build("PyCall")
	using PyCall
	
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,  padding=(0, 10, 0, 0)),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS = Makie.wong_colors()
end

# ╔═╡ 0132b4ab-0447-4546-b412-ec598b20d21d
md"""
# Retrievals

[Intro text here]

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ db524678-9ee2-4934-b1bb-6a2f13bf0fa6
dirpath = 
"data/retrievals/" *
"spot_lower_fit_R0/"*
"WASP50_E1_NoHet_FitP0_NoClouds_NoHaze_fitR0_Na_K"

# ╔═╡ 38f37547-9663-464d-9983-4fd3bdbc79b6
retr = CSV.File("$(dirpath)/retr_Magellan_IMACS.txt")

# ╔═╡ b245f244-897e-4eb2-a6a6-cef7c85ca390
retr_model = CSV.File("$(dirpath)/retr_model.txt") 

# ╔═╡ 66ad752e-29a2-4d1d-8985-4bda58138e31
retr_model_sampled = CSV.File("$(dirpath)/retr_model_sampled_Magellan_IMACS.txt")

# ╔═╡ e801501c-a882-4f2d-bbc1-40028c1c91d8
begin
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	errorbars!(ax, retr.wav, retr.flux, retr.flux_err)
	scatter!(ax, retr.wav, retr.flux)
	
	lines!(ax, retr_model.wav, retr_model.flux)
	scatter!(ax, retr_model_sampled.wav, retr_model_sampled.flux)
	
	xlims!(ax, 0.2, 1.1)
	
	fig
end

# ╔═╡ 1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
md"""
## Notebook setup
"""

# ╔═╡ Cell order:
# ╟─0132b4ab-0447-4546-b412-ec598b20d21d
# ╠═db524678-9ee2-4934-b1bb-6a2f13bf0fa6
# ╠═38f37547-9663-464d-9983-4fd3bdbc79b6
# ╠═b245f244-897e-4eb2-a6a6-cef7c85ca390
# ╠═66ad752e-29a2-4d1d-8985-4bda58138e31
# ╠═e801501c-a882-4f2d-bbc1-40028c1c91d8
# ╟─1eff1230-2423-4ac3-8e9b-f4e7bcd0121b
# ╠═239a91a6-f68a-11eb-14fd-0ba8d08b08f9
