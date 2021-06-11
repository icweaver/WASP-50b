### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 89e9971c-7e84-43d5-95c7-b9cfc3b1b91f
using TimeSeries

# ╔═╡ 5de1e52e-8074-4d38-9143-f8f3cd14f4c5
using Plots

# ╔═╡ 75624bac-af3f-47d1-9462-7840bac427a8
using Dates

# ╔═╡ 9e2ce576-c9bd-11eb-0699-47af13e79589
begin
	using AlgebraOfGraphics
	using CSV
	#using CairoMakie
	using Colors
	using DataFrames
	using DataFramesMeta
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
# Photometric monitoring

$(TableOfContents(title="📖 Table of Contents"))
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
$(join(name.(extrema(df_phot_mon.hjd)), " - ")) from two cameras (bd: 2013 December, bh: 2015 July), both in the V band. Since we will be using both cameras, let's sort the entire dataframe by time and plot the results next.
"""

# ╔═╡ 92548af3-9a26-4202-88f2-ba3a31181686
begin
	df_sorted = sort(df_phot_mon, :hjd)
	t, f, f_err = eachcol(df_sorted[!, [:hjd, :mag, :mag_err]])
end

# ╔═╡ 1f62da17-8e2b-4420-ab50-7d28c80ebfed
plotly()

# ╔═╡ 886fe663-ea95-47e1-97d2-7a467ab41394
dates = Date(2018, 1, 1):Day(1):Date(2018, 12, 31)

# ╔═╡ 2d43eb74-1f56-4a0f-bd96-865cb9dc35c5
ta = TimeArray(julian2datetime.(t), f)

# ╔═╡ 8ee19677-e72b-441d-9429-687ea246cf60
begin
	p = scatter(ta)
	scatter!(p, moving(mean, ta, 13))
end

# ╔═╡ 9f99c140-a1ba-488c-8736-af30db45cc02
moving(mean, ta, 5) |> length

# ╔═╡ 3ae4c79f-649e-4055-932f-a8870a335134
julian2datetime.(t)

# ╔═╡ a8b37aba-7613-4019-a051-6fda0f8d8686
plot

# ╔═╡ 53ec6a2a-1f7e-44f8-b5c5-2f90a4960542
diff(t) |> mean

# ╔═╡ 726571bd-0d4b-453e-92f5-a5a123974f02
t_roll = rolling(mean, t, 2)

# ╔═╡ d94e6e76-3e78-45f6-9681-0f5b345b0831
f_roll = rolling(mean, f, 2)

# ╔═╡ 278701a9-5503-4bde-9416-af8e84f13923
@which mean

# ╔═╡ be05a66d-bc4f-42ef-9734-98c9a8ff2b88
let
	fig = figure
end

# ╔═╡ 2eb8be9c-bf51-45c3-b44d-3028344bd13e
utc_transit_dates = ["2013-12-19", "2015-09-27", "2016-12-11"]

# ╔═╡ 94506bca-159c-4783-82c5-9837527f1c6e
transit_dates = DateTime.(utc_transit_dates)

# ╔═╡ 7dd2f68b-8f30-405d-bfa6-ca6542dd64b1
julian_transit_dates = datetime2julian.(transit_dates) |> collect

# ╔═╡ 8e6008ca-a762-4450-a09e-bd0c2bbac4f2
begin
	set_aog_theme!()
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
		
	fig |> as_svg
end

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

# ╔═╡ 225deb3e-71a5-4fd2-b49e-64ae55c1bf05
let
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="Time (HJD)", ylabel="Magnitude")
	
	errorbars!(ax, t, f, f_err, color=COLORS[1])
	scatter!(ax, t, f)
	
	scatter!(ax, t_roll, f_roll)
	
	fig
end

# ╔═╡ ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╟─670b88e4-1e96-494d-bfcc-235092bb6e96
# ╠═fa233a2c-7e89-4e71-85b8-824c5c341650
# ╠═9094c4a4-3b75-4e21-97a7-600de734867b
# ╠═2b2dd83b-ce99-4551-b8aa-79ca6db5dd06
# ╟─6eaf882c-0cb5-415f-b8fe-c071ee25a895
# ╠═92548af3-9a26-4202-88f2-ba3a31181686
# ╠═225deb3e-71a5-4fd2-b49e-64ae55c1bf05
# ╠═89e9971c-7e84-43d5-95c7-b9cfc3b1b91f
# ╠═1f62da17-8e2b-4420-ab50-7d28c80ebfed
# ╠═8ee19677-e72b-441d-9429-687ea246cf60
# ╠═9f99c140-a1ba-488c-8736-af30db45cc02
# ╠═886fe663-ea95-47e1-97d2-7a467ab41394
# ╠═2d43eb74-1f56-4a0f-bd96-865cb9dc35c5
# ╠═3ae4c79f-649e-4055-932f-a8870a335134
# ╠═a8b37aba-7613-4019-a051-6fda0f8d8686
# ╠═5de1e52e-8074-4d38-9143-f8f3cd14f4c5
# ╠═53ec6a2a-1f7e-44f8-b5c5-2f90a4960542
# ╠═726571bd-0d4b-453e-92f5-a5a123974f02
# ╠═d94e6e76-3e78-45f6-9681-0f5b345b0831
# ╠═278701a9-5503-4bde-9416-af8e84f13923
# ╠═8e6008ca-a762-4450-a09e-bd0c2bbac4f2
# ╠═be05a66d-bc4f-42ef-9734-98c9a8ff2b88
# ╠═2eb8be9c-bf51-45c3-b44d-3028344bd13e
# ╠═7dd2f68b-8f30-405d-bfa6-ca6542dd64b1
# ╠═94506bca-159c-4783-82c5-9837527f1c6e
# ╠═75624bac-af3f-47d1-9462-7840bac427a8
# ╟─7ffcd329-91ce-4936-b86c-b9c11aa07a2b
# ╠═3847c155-0ef4-463b-abfd-5724cf0d28bb
# ╟─ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
# ╠═9e2ce576-c9bd-11eb-0699-47af13e79589
