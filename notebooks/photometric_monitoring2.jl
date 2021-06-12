### A Pluto.jl notebook ###
# v0.14.8

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

# ╔═╡ fbed4840-cada-11eb-0c6c-d38d9654159f
begin
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using Colors
	using DataFrames
	using DataFramesMeta
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using KernelDensity
	using Latexify
	using Measurements
	using NaturalSort
	using OrderedCollections
	using Printf
	using PyCall
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
end

# ╔═╡ 6fd6be4c-0de1-4f12-8cd3-b773eb3666f2
md"""
# Photometric Monitoring 2 ⭐

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ d8f23d33-3ac3-417b-a2f5-cf7dc61beb01
md"""
## Load data
"""

# ╔═╡ 3cec41d5-bb87-446b-b39e-80b949d108da
md"""
Let's start by loading in the data and summarizing its contents:
"""

# ╔═╡ 4834702f-7dd6-4e9b-bfb7-8c9c2ba7651b
df_phot_mon = CSV.File(
	"data/photometric/AP37847073.csv",
	normalizenames = true,
) |> DataFrame

# ╔═╡ ed9ea9e6-9728-4313-b082-7d8515a24591
describe(df_phot_mon, :all)

# ╔═╡ 771bb36b-affa-491b-b732-2bfabd76d8ef
function name(dt_hjd)
	dt = julian2datetime(dt_hjd)
	return "$(year(dt)) $(monthname(dt)) $(day(dt))"
end

# ╔═╡ 04761ee0-2f91-4b4c-9100-412f150d3cf2
md"""
According to the table above, the data spans from
$(join(name.(extrema(df_phot_mon.hjd)), " - ")) from two cameras (bd: 2013 December, bh: 2015 July), both in the V band. Since we will be using both cameras, let's sort the entire dataframe by time and plot the results next.
"""

# ╔═╡ 881a2bbf-cb5b-4242-af71-19e876720dec
begin
	df_sorted = sort(df_phot_mon, :hjd)
	t, f, f_err = eachcol(df_sorted[!, [:hjd, :mag, :mag_err]])
end

# ╔═╡ 533d8696-425a-4788-95ef-1370d548ada3
md"""
## Plot
"""

# ╔═╡ 1d789b40-8fa1-4983-aac2-7275e5a0e29c
md"""
We show the original and averaged data below:
"""

# ╔═╡ 7cb85f49-b3c6-4145-add2-6719e2c7201e
md"""
Δt = $(@bind Δt Slider(1:10, default=7, show_value=true)) days
"""

# ╔═╡ 71de9258-8b02-4a88-951a-846736cd4407
md"""
## Plot description
"""

# ╔═╡ 6db1d39c-c03c-461b-9f1e-4353948dd83a
md"""
We produced the above plot by performing a rolling average over sub-groups of the data. Each sub-group contains timeseries points within $Δt$ of the first point in that sub-group. We then move on to the next group and perform the same grouping operation:

!!! note
	TODO: Find a good Δt to use
"""

# ╔═╡ 81a40d75-0880-47fa-b2c6-5b79c887d227
# Given a vector `v`, return the set of groups satisfying the condition that all points in that group are ≤ Δv of the first point
function group_vector(v, Δv)
	groups = Vector{Float64}[]
	current_idx = 0
	group_idxs = Vector{Int}[]
	v_idxs = 1:length(v)
	while current_idx < length(v)
		is_group = @. 0.0 ≤ v - v[begin+current_idx] ≤ Δv
		group = v[is_group]
		group_idx = v_idxs[is_group]
		push!(groups, group)
		push!(group_idxs, group_idx)
		current_idx += length(group)
	end
	return groups, group_idxs
end

# ╔═╡ 3a53aee7-a4fb-4275-aa80-9d836f84c795
md"""
Once the grouping was complete, we then averaged over each group to produce the filtered points above:
"""

# ╔═╡ decb01e5-16d7-4c88-a202-5ad93acfbcdf
# Averages the xs and associated ys in each grouping of x
function bin_data(x, y; Δ=0.1)
	x_groups, group_idxs = group_vector(x, Δ)
	y_groups = [y[idxs] for idxs in group_idxs]
	x_binned = mean.(x_group for x_group ∈ x_groups)
	y_binned = mean.(y_group for y_group ∈ y_groups)
	return x_binned, y_binned
end

# ╔═╡ 7b11cbfa-1f56-41cf-b622-fdf779a61a98
t_binned, f_binned = bin_data(t, f, Δ=Δt);

# ╔═╡ 0fba0bbf-1290-47ef-904a-5b576b70d07c
md"""
## Plot configs
"""

# ╔═╡ 02277a5e-efba-48a5-996e-ba639574e3f2
begin
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	#const COLORS = to_colormap(:seaborn_colorblind6, 8)[[8, 6, 4, 1]]
	const COLORS = parse.(Colorant,
		[
			"#a6cee3",  # Cyan
			"#fdbf6f",  # Yellow
			"#ff7f00",  # Orange
			"#1f78b4",  # Blue
			"slategray",
			# "plum",
			# "#956cb4",  # Purple
			# "mediumaquamarine",
			# "#029e73",  # Green
		]
	)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			palette = (color=COLORS,),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	
	COLORS
end

# ╔═╡ 81ff1d9b-f4e4-4c18-8570-9584f04cecef
let
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="Time (HJD)", ylabel="Magnitude")
	
	scatter!(ax, t, f, markersize=15, color=(COLORS[end], 0.5))
	scatter!(ax, t_binned, f_binned, color=COLORS[3])
	
	fig
end

# ╔═╡ cba385b1-7d69-408a-a622-97e6498ae978
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╟─6fd6be4c-0de1-4f12-8cd3-b773eb3666f2
# ╟─d8f23d33-3ac3-417b-a2f5-cf7dc61beb01
# ╟─3cec41d5-bb87-446b-b39e-80b949d108da
# ╠═4834702f-7dd6-4e9b-bfb7-8c9c2ba7651b
# ╠═ed9ea9e6-9728-4313-b082-7d8515a24591
# ╟─04761ee0-2f91-4b4c-9100-412f150d3cf2
# ╟─771bb36b-affa-491b-b732-2bfabd76d8ef
# ╠═881a2bbf-cb5b-4242-af71-19e876720dec
# ╟─533d8696-425a-4788-95ef-1370d548ada3
# ╟─1d789b40-8fa1-4983-aac2-7275e5a0e29c
# ╟─7cb85f49-b3c6-4145-add2-6719e2c7201e
# ╠═81ff1d9b-f4e4-4c18-8570-9584f04cecef
# ╠═7b11cbfa-1f56-41cf-b622-fdf779a61a98
# ╟─71de9258-8b02-4a88-951a-846736cd4407
# ╟─6db1d39c-c03c-461b-9f1e-4353948dd83a
# ╠═81a40d75-0880-47fa-b2c6-5b79c887d227
# ╟─3a53aee7-a4fb-4275-aa80-9d836f84c795
# ╠═decb01e5-16d7-4c88-a202-5ad93acfbcdf
# ╟─0fba0bbf-1290-47ef-904a-5b576b70d07c
# ╠═02277a5e-efba-48a5-996e-ba639574e3f2
# ╟─cba385b1-7d69-408a-a622-97e6498ae978
# ╠═fbed4840-cada-11eb-0c6c-d38d9654159f
