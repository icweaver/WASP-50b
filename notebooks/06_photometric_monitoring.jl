### A Pluto.jl notebook ###
# v0.15.0

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

# ╔═╡ 9e2ce576-c9bd-11eb-0699-47af13e79589
begin
	import Pkg
	Pkg.activate(joinpath(@__DIR__, ".."))
	
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

# ╔═╡ 670b88e4-1e96-494d-bfcc-235092bb6e96
md"""
# Photometric Monitoring

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 0cbe4263-799f-4ee3-9a94-3ba879528b01
md"""
## ASAS-SN
"""

# ╔═╡ fa233a2c-7e89-4e71-85b8-824c5c341650
df_phot_mon = CSV.File(
	"data/photometric/AP37847073.csv",
	normalizenames = true,
) |> DataFrame

# ╔═╡ 9094c4a4-3b75-4e21-97a7-600de734867b
describe(df_phot_mon, :all)

# ╔═╡ 2b2dd83b-ce99-4551-b8aa-79ca6db5dd06
function name(dt_hjd)
	dt = julian2datetime(dt_hjd)
	return "$(year(dt)) $(monthname(dt)) $(day(dt))"
end

# ╔═╡ 6eaf882c-0cb5-415f-b8fe-c071ee25a895
md"""
According to the table above, the data spans from
$(join(name.(extrema(df_phot_mon.hjd)), " - ")) from two cameras (bd: 2013 December, bh: 2015 July), both in the V band. Since we will be using both cameras, let's sort the entire dataframe by time and plot the results next:
"""

# ╔═╡ 92548af3-9a26-4202-88f2-ba3a31181686
begin
	df_sorted = sort(df_phot_mon, :hjd)
	t_ASASSN, f_ASASSN, f_err_ASASSN = eachcol(df_sorted[!, [:hjd, :mag, :mag_err]])
	
	utc_transit_dates = ["2013-12-19", "2015-09-27", "2016-12-11"]
	transit_dates = DateTime.(utc_transit_dates)
	julian_transit_dates = datetime2julian.(transit_dates) |> collect
end;

# ╔═╡ 8e6008ca-a762-4450-a09e-bd0c2bbac4f2
let
	phot_mon = data(df_phot_mon) * mapping(
	    :hjd => "Time (HJD)",
	    :mag => "Magnitude",
		color = :camera,
	)
	
	fig = draw(phot_mon)
	ax = current_axis()
	
	vlines!(ax, julian_transit_dates;
		linestyle = :dash,
		color = :darkgrey,
	)
	
	for (i, (utc, jd)) in enumerate(zip(utc_transit_dates, julian_transit_dates))
		text!(ax, "Transit $i\n$utc";
			position = Point2f0(jd, 11.8),
			textsize = 14,
			align = (:left, :center),
			offset = Point2f0(10, 0),
			color = :grey,
		)
	end
		
	#save("phot_mon.png", fig, px_per_unit=3)
	
	fig #|> as_svg
end

# ╔═╡ dbe317fe-540d-44e8-b8e7-6c465c79559f
md"""
Δt = $(@bind Δt Slider(1:10, default=7, show_value=true)) days
"""

# ╔═╡ 78d85c7c-da06-41ab-915b-48d93a010967
md"""
## TESS
"""

# ╔═╡ e7eba854-3846-4ba4-bbdb-c5f7dbdf08a7
df_TESS_S04 = CSV.File(
	"data/photometric/TESS_baseline_sector_04.csv",
) |> DataFrame

# ╔═╡ e1648b37-b52b-4c46-a5fb-4e5bc2743cd4
df_TESS_S31 = CSV.File(
	"data/photometric/TESS_baseline_sector_31.csv",
) |> DataFrame

# ╔═╡ 09c666f7-a8b4-4b47-b9fd-1351e8bd28f9
Δ_TESS = 200.0 * (1/86_400) * 10

# ╔═╡ 5cca08d9-e035-4a15-b575-069e6b89b6db
t_TESS_S04, f_TESS_S04, f_err_TESS_S04 = eachcol(
	df_TESS_S04[!, [:time, :flux, :flux_err]]
)

# ╔═╡ a87d6163-4fd9-49c3-a83b-ab9c727edf99
t_TESS_S31, f_TESS_S31, f_err_TESS_S31 = eachcol(
	df_TESS_S31[!, [:time, :flux, :flux_err]]
)

# ╔═╡ 18223d42-66d8-40d1-9d89-be8af46853e2
md"""
## Helper Functions
"""

# ╔═╡ 976bd30f-2bae-4c5d-93dc-70b4df303840
md"""
We produced the above plots by performing a rolling average over sub-groups of the data. Each sub-group contains timeseries points within $Δt$ of the first point in that sub-group. We then move on to the next group and perform the same grouping operation:

!!! note
	TODO: Find a good Δt to use
"""

# ╔═╡ 28ce86fd-9679-4775-80ff-75479baa8a0c
# Given a sorted vector `v`, return the set of groups (by idx) satisfying the condition that all points in that group are ≤ Δv of the first point
function idxs_of_groups(v, Δv)
	idx_start = 1
	groups = Vector{Int}[]
	group = [idx_start]
	for i in 2:length(v)
		if v[i] - v[idx_start] ≤ Δv
			push!(group, i)
		else
			push!(groups, group)
			idx_start = i
			group = [i]
		end
	end
	push!(groups, group)
	
	return groups
end

# ╔═╡ ed0feb32-85a2-4d3b-bb96-47eb18e29dd3
# Averages the xs and associated ys in each grouping of x
function bin_data(x, y; Δ=0.1)
	idxs = idxs_of_groups(x, Δ)
	groups = [(x[idx], y[idx]) for idx in idxs]
	binned = [(mean(group[1]), weightedmean(group[2]).val) for group ∈ groups]
	return binned
end

# ╔═╡ 066837a0-6390-4832-85b8-2143f7ddb471
binned_vals_ASAS_SN = bin_data(t_ASASSN, f_ASASSN .± f_err_ASASSN, Δ=Δt);

# ╔═╡ 398f0134-f6be-4d4b-be99-1bfc03e81d66
binned_vals_TESS_S04 = bin_data(
	t_TESS_S04, f_TESS_S04 .± f_err_TESS_S04, Δ=Δ_TESS
)

# ╔═╡ 4a434497-bfe4-4033-9bbe-9b36783705ea
binned_vals_TESS_S31 = bin_data(
	t_TESS_S31, f_TESS_S31 .± f_err_TESS_S31, Δ=Δ_TESS
)

# ╔═╡ 7ffcd329-91ce-4936-b86c-b9c11aa07a2b
md"""
## Plot configs
"""

# ╔═╡ 3847c155-0ef4-463b-abfd-5724cf0d28bb
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

# ╔═╡ 7370a1d9-4f8e-4788-adac-b8be2bcc9643
function plot_phot(ax, t, f, f_err, binned_vals)
	# Original data
	errorbars!(ax, t, f, f_err, color=(COLORS[end], 0.25))
	scatter!(ax, t, f, color=(COLORS[end], 0.25))
	
	# Averaged data
	scatter!(ax, binned_vals, markersize=14, color=COLORS[end-1])
end

# ╔═╡ 2ee9b0ca-f6a2-47db-abfc-9a44e64b2e42
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	vlines!(ax, julian_transit_dates;
		linestyle = :dash,
		color = :darkgrey,
	)
	
	for (i, (utc, jd)) in enumerate(zip(utc_transit_dates, julian_transit_dates))
		text!(ax, "Transit $i\n$utc";
			position = Point2f0(jd, 11.8),
			textsize = 14,
			align = (:left, :center),
			offset = Point2f0(10, 0),
			color = :grey,
		)
	end
	
	plot_phot(ax, t_ASASSN, f_ASASSN, f_err_ASASSN, binned_vals_ASAS_SN)
	fig
end

# ╔═╡ 8906a2a2-65c9-4dc1-aaef-078a6ddaaff2
let
	fig = Figure()
	ax_top = Axis(fig[1, 1],)
	ax_bottom = Axis(fig[2, 1])

	plot_phot(ax_top, t_TESS_S04, f_TESS_S04, f_err_TESS_S04, binned_vals_TESS_S04)
	plot_phot(ax_bottom, t_TESS_S31, f_TESS_S31, f_err_TESS_S31, binned_vals_TESS_S31)
	
	fig #|> as_svg
end

# ╔═╡ ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
md"""
## Packages
"""

# ╔═╡ f8e37bb8-bdd8-4eea-82c3-1b1f3a5617a1
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
# ╟─670b88e4-1e96-494d-bfcc-235092bb6e96
# ╟─0cbe4263-799f-4ee3-9a94-3ba879528b01
# ╠═fa233a2c-7e89-4e71-85b8-824c5c341650
# ╠═9094c4a4-3b75-4e21-97a7-600de734867b
# ╠═2b2dd83b-ce99-4551-b8aa-79ca6db5dd06
# ╟─6eaf882c-0cb5-415f-b8fe-c071ee25a895
# ╠═92548af3-9a26-4202-88f2-ba3a31181686
# ╠═8e6008ca-a762-4450-a09e-bd0c2bbac4f2
# ╟─dbe317fe-540d-44e8-b8e7-6c465c79559f
# ╠═2ee9b0ca-f6a2-47db-abfc-9a44e64b2e42
# ╠═066837a0-6390-4832-85b8-2143f7ddb471
# ╟─78d85c7c-da06-41ab-915b-48d93a010967
# ╠═e7eba854-3846-4ba4-bbdb-c5f7dbdf08a7
# ╠═e1648b37-b52b-4c46-a5fb-4e5bc2743cd4
# ╠═09c666f7-a8b4-4b47-b9fd-1351e8bd28f9
# ╠═5cca08d9-e035-4a15-b575-069e6b89b6db
# ╠═398f0134-f6be-4d4b-be99-1bfc03e81d66
# ╠═a87d6163-4fd9-49c3-a83b-ab9c727edf99
# ╠═4a434497-bfe4-4033-9bbe-9b36783705ea
# ╠═8906a2a2-65c9-4dc1-aaef-078a6ddaaff2
# ╟─18223d42-66d8-40d1-9d89-be8af46853e2
# ╟─976bd30f-2bae-4c5d-93dc-70b4df303840
# ╠═28ce86fd-9679-4775-80ff-75479baa8a0c
# ╠═ed0feb32-85a2-4d3b-bb96-47eb18e29dd3
# ╠═7370a1d9-4f8e-4788-adac-b8be2bcc9643
# ╟─7ffcd329-91ce-4936-b86c-b9c11aa07a2b
# ╠═3847c155-0ef4-463b-abfd-5724cf0d28bb
# ╟─ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
# ╠═9e2ce576-c9bd-11eb-0699-47af13e79589
# ╟─f8e37bb8-bdd8-4eea-82c3-1b1f3a5617a1
