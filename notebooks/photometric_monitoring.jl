### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ 75624bac-af3f-47d1-9462-7840bac427a8
using Dates

# ╔═╡ 9e2ce576-c9bd-11eb-0699-47af13e79589
begin
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
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
From the table above, we see that the data spans from
$(join(name.(extrema(df_phot_mon.hjd)), " - ")) from two cameras (bd: 2013 December, bh: 2015 July), both in the V band:
"""

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
# ╠═8e6008ca-a762-4450-a09e-bd0c2bbac4f2
# ╠═2eb8be9c-bf51-45c3-b44d-3028344bd13e
# ╠═7dd2f68b-8f30-405d-bfa6-ca6542dd64b1
# ╠═94506bca-159c-4783-82c5-9837527f1c6e
# ╠═75624bac-af3f-47d1-9462-7840bac427a8
# ╟─7ffcd329-91ce-4936-b86c-b9c11aa07a2b
# ╠═3847c155-0ef4-463b-abfd-5724cf0d28bb
# ╟─ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
# ╠═9e2ce576-c9bd-11eb-0699-47af13e79589
