### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# ╔═╡ 24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics, CairoMakie, CSV, DataFrames, Unitful, UnitfulAstro
	using DataFrameMacros
	using HTTP
	using Unitful
	using PlutoUI
	using Unitful: k, G
	using Chain
	using NaturalSort
	
	set_aog_theme!()
	
	update_theme!(
		Theme(
			Axis = (xlabelsize=18, ylabelsize=18,),
			Label = (textsize=18,  padding=(0, 10, 0, 0)),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Text = (font=AlgebraOfGraphics.firasans("Light"),),
			fontsize = 18,
			rowgap = 0,
			colgap = 0,
		)
	)
	COLORS = Makie.wong_colors()
end

# ╔═╡ cd13d7f3-0ea3-4631-afd9-5f3e359000e6
md"""
# HGJH and CO/CH₄ population

In this notebook we will explore the possible targets amenable to atmopsheric characterization near:
1. the high-gravity hot-Jupiter (HGHJ) parameter space (Tₚ ∼ 1000 K, g ∼ 30 m/s²)
1. the potential molecular transition between CO – CH₄ near the 1200 K boundary predicted from analagous studies of brown dwarfs (BDs)

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 6b06701b-05e2-4284-a308-e9edeb65a648
md"""
## Data sample

We will draw our sample of each sample above from the archive of transiting exoplanet data from ground and space-based observations, which we query from the NASA Exoplanet Archive [TAP API](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html):
"""

# ╔═╡ f396cda3-f535-4ad9-b771-7ccbd45c54f3
df_all = let
		columns = [
		"pl_name",
		"disc_facility",
		"tic_id",
		"pl_rade",
		"pl_bmasse",
		"pl_orbsmax",
		"st_teff",
		"st_rad",
		"sy_jmag",
	]
	url = "https://exoplanetarchive.ipac.caltech.edu/TAP"
	#cond = "tran_flag+=1+and+pl_eqt+<+1000+and+pl_rade+<+4"
	cond = "tran_flag+=1"
	query = "select+$(join(columns, ','))+from+pscomppars+where+$(cond)&format=csv"
	request = HTTP.get("$(url)/sync?query=$(query)")
	CSV.read(request.body, DataFrame)
end

# ╔═╡ 4d1a7740-24c7-4cec-b788-a386bc25f836
md"""
We next compute some relevant quantities from this table to help organize each population:
"""

# ╔═╡ e0365154-d6c8-4db2-bb85-bf2536a3aa74
function compute_Teq(T, R, a; α)
	T * (1.0 - α)^0.25 * sqrt(0.5 * R / a)
end

# ╔═╡ 9aed232f-ec74-4ec6-9ae7-06b90539833b
# begin
# 	df = @chain df_all begin
# 		dropmissing([:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag, :pl_trandep])
# 		@aside g_SI = let
# 			g = compute_g.(_.pl_bmasse*u"Mearth", _.pl_rade*u"Rearth")
# 			ustrip.(u"m/s^2", g)
# 		end
# 		@aside H_km = let
# 			H = compute_H.(2.0*u"u", _.pl_eqt*u"K", g_SI*u"m/s^2")
# 			ustrip.(u"km", H)
# 		end
# 		@transform begin
# 			:g_SI = g_SI
# 			:H_km = H_km
# 			:TSM = compute_TSM.(
# 				:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag
# 			)
# 		end
# 		@transform :TSMR = :TSM ./ :TSM[:pl_name .== "HAT-P-23 b"][1]
# 	end
# end

# ╔═╡ 4cbbb1e8-e5fb-4ab0-a7e6-7881c2dde032
md"""
With this list of $(nrow(df_all)) transiting exoplanets with observed stellar and planetary parameters, we now turn to visualizing the subset of each population from the list above.
"""

# ╔═╡ 7f956b36-ce65-4e4e-afa5-9b97b9e06954
md"""
## HGHJ population

We define this population to be the sample of transiting exoplanets with ... We defined our relative TSM (TSMR) to be based on the TSM of HAT-P-23b for comparison.
"""

# ╔═╡ 0f9262ef-b774-45bc-bdab-46860779683d
md"""
!!! note
	Inspired from [warm-worlds](https://github.com/nespinoza/warm-worlds)
"""

# ╔═╡ 18094afc-b77f-4cae-a30c-2691d34125d8
md"""
!!! warning "TODO"
* WASP-33: active star, difficult analysis
* TOI 1581: double check
* WASP-4 : Already in ACCESS survey, can re-analyze
* Focus on top 6-12 targets, review targets in literature, what can be done from the South?
"""

# ╔═╡ 958453c3-7993-4620-ab7f-e7ad79781dd5
val(df, name, col) = df[df.pl_name .== name, col][1]

# ╔═╡ f07ad06b-81d2-454f-988f-a7ae1713eac4
function annotate_text!(ax, t, p1, p2, l, lm; align=(:center, :center))
	hyp = 2.0*lm + l
	eps = lm*(p2 - p1) / hyp
	text!(ax, t, position=p1, align=align)
	lines!(ax, [p1 + eps, p2 - eps], color=:darkgrey, linewidth=1)
end

# ╔═╡ 2776646e-47e7-4b9e-ab91-4035bc6df99f
function compute_scale_factor(Rp)
	if Rp < 1.5
		f = 0.190
	elseif 1.5 ≤ Rp 2.75
		f = 1.26
	elseif 2.75 ≤ Rp < 4.0
		f = 1.28
	elseif 4.0 ≤ Rp < 10.0
		f = 1.15
	else
		f = 1.15 # Only goes up to 1O R_earth in Kempton+2018
	end
	return f
end

# ╔═╡ c7960066-cc33-480c-807b-c56ead4262bf
# compute TSM, assuming M, R in Earth units, T in Kelvin
function compute_TSM(Rp, Teq, Mp, Rs, J; denom=1.0)
	f = compute_scale_factor(Rp)
	return f * (Rp^3 * Teq / (Mp * Rs^2)) * 10.0^(-J/5.0) / denom
end

# ╔═╡ abaee9cc-9841-4b6b-ad33-2093c27422c8
compute_g(M, R) = G * M / R^2

# ╔═╡ 1e587a84-ed43-4fac-81cf-134a4f3d65d5
compute_H(μ, T, g) = k * T / (μ * g)

# ╔═╡ c1e63cf3-7f30-4858-bdd6-125d2a99529f
compute_ΔD(H, Rₚ, Rₛ) = 2.0 * H * Rₚ/Rₛ^2

# ╔═╡ 2702ba80-993c-4cca-bd14-0d4d6b67362b
begin
	df = dropmissing(
		df_all, [:pl_rade, :pl_bmasse, :pl_orbsmax, :st_teff, :sy_jmag]
	)
	df.pl_eqt = let
		Teq = compute_Teq.(df.st_teff*u"K", df.st_rad*u"Rsun", df.pl_orbsmax*u"AU"; α=0.1)
		ustrip.(u"K", Teq)
	end
	df.g_SI = let
		g = compute_g.(df.pl_bmasse*u"Mearth", df.pl_rade*u"Rearth")
		ustrip.(u"m/s^2", g)
	end
	df.H_km = let
		H = compute_H.(2.0*u"u", df.pl_eqt*u"K", df.g_SI*u"m/s^2")
		ustrip.(u"km", H)
	end
	df.TSM = compute_TSM.(df.pl_rade, df.pl_eqt, df.pl_bmasse, df.st_rad, df.sy_jmag)
	df.TSMR = df.TSM ./ df.TSM[df.pl_name .== "HAT-P-23 b"][1]
	df.ΔD_ppm = let
		ΔD = compute_ΔD.(df.H_km*u"km", df.pl_rade*u"Rearth", df.st_rad*u"Rsun")
		uconvert.(NoUnits, ΔD) * 5.0 * 1e6
	end
	df.ΔDR_ppm = df.ΔD_ppm ./ df.ΔD_ppm[df.pl_name .== "HAT-P-23 b"][1]
	df
end

# ╔═╡ c98c5618-4bac-4322-b4c3-c76181f47889
df_HGHJs_all = @chain df begin
	@subset (1.0 ≤ :TSMR) &
	(15.0 ≤ :g_SI ≤ 53) &
	(:pl_eqt ≤ 4900.0)
	sort(:TSMR, rev=true)
end

# ╔═╡ d62b5506-1411-49f2-afe3-d4aec70641a1
df_HGHJs = @subset(
	df_HGHJs_all, :pl_name .∈ Ref(["HAT-P-23 b", "WASP-43 b", "WASP-50 b"])
)

# ╔═╡ c1cd9292-28b9-4206-b128-608aaf30ff9c
# TODO: Place latitude constraints
let
	fig = Figure(resolution=(800, 700))
	ax = Axis(fig[1, 1], limits=((0, nothing), (10, 55)))

	# HGHJ g boundary
	hlines!(ax, 20.0, color=:darkgrey, linestyle=:dash)
	
	# Phase plot
	markersize_factor = 10.0
	m = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:g_SI => "Surface gravity (m/s²)",
		color = :ΔD_ppm => "ΔD (ppm)",
		markersize = :TSMR => (x -> markersize_factor*x),
	)
	m2 = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:g_SI => "Surface gravity (m/s²)",
		#color = :ΔD_ppm => "ΔD (ppm)",
		#markersize = :TSMR => (x -> markersize_factor*x),
	)
	marker_open = visual(marker='○') # ○ ●
	marker_closed = visual(markersize=15, marker='+', color=:white)
	plt = m*data(df_HGHJs_all)
	fg = draw!(ax, plt)
	colorbar!(fig[1, 2], fg)
	
	# HGHJ g boundary
	# ax = fg.grid[1, 1].axis
	#hlines!(ax, 20.0, color=:darkgrey, linestyle=:dash)
	
	# Annotate HGHJs with tspec observations
	HP23x, HP23y = val.(Ref(df_HGHJs), Ref("HAT-P-23 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"HAT-P-23 b",
		Point2f((HP23x[1], HP23y[1])) .- (-600, -3),
		Point2f((HP23x[1], HP23y[1])),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	W43x, W43y = val.(Ref(df_HGHJs), Ref("WASP-43 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"WASP-43 b",
		Point2f((W43x[1], W43y[1])) .- (800, 5),
		Point2f((W43x[1], W43y[1])),
		0.3,
		0.1;
		align = (:center, :top),
	)
	W50x, W50y = val.(Ref(df_HGHJs), Ref("WASP-50 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"WASP-50 b",
		Point2f((W50x[1], W50y[1])) .- (-300, -8),
		Point2f((W50x[1], W50y[1])),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	
	# TSMR legend
	tsmrs = [14, 4, 1]
	axislegend(
		ax,
		[MarkerElement(marker='○', markersize=markersize_factor*ms) for ms ∈ tsmrs],
		["$tsmr" for tsmr ∈ tsmrs],
		"TSMR",
		position = :rt,
		patchsize = (120, 100),
		framevisible = true,
		padding = (5, 5, -24, 10),
		margin = (10, 10, 0, 0),
		titlegap = 24,
	)

	path = "../../ACCESS_WASP-50b/figures/pop"
	mkpath(path)
	save("$(path)/hg_pop.png", fig)
	
	fig
end

# ╔═╡ d2d6452b-426c-4c61-8dc4-d210c0cd9d41
@which println(1000000.)

# ╔═╡ 2476f26e-f7cc-4f8e-ac66-60b85a46cb2c
md"""
## CO – CH₄ transition
"""

# ╔═╡ 8e4b1149-abf1-4ded-b20f-4765d2ee49a9
md"""
For this population, we will make a similar plot and compare to the current list of ACCESS targets that straddle this boundary.
"""

# ╔═╡ 8cece13d-34cb-40df-8986-ac0dad210e58
df_ACCESS = innerjoin(CSV.read("data/pop/ACCESS.csv", DataFrame), df, on=:pl_name)

# ╔═╡ df499885-3f7d-4293-bcf5-9c761ac3ccd2
sort(df_ACCESS, :pl_name, lt=natural)

# ╔═╡ d6791747-7503-49a9-903c-c479fc0c3d49
species = (
	"H₂O" => (273.15, (:center, :center)),
	"NH₃" => (583.0, (:right, :center)),
	"N₂" => (583.0, (:left, :baseline)),
	"CH₄" => (1000.0, (:right, :center)),
	"CO" => (1000.0, (:left, :baseline)),
	"MnS" => (1350.0, (:center, :center)),
)

# ╔═╡ 60016c3f-6968-4d3c-ac73-d324b2a071e0
begin
	fig = Figure()
	ax = Axis(fig[1, 1], limits=(0, 3_600, 0, nothing))
	
	# Condensation temps
	for (i, (name, (T, align))) in enumerate(species)
		vlines!(ax, T, color=:darkgrey, linewidth=1.0, linestyle=:dash)
		text!(ax, name, position=(T, 16.0 + i), align=align, textsize=16)
	end
	vspan!(ax, 1625.0, 1875.0, color=(:darkgrey, 0.25))
	text!(ax, "Silicates/Metal-oxides";
		position = (0.5*(1875.0+1625.0), 22.0),
		textsize = 16,
	)
	
	# Phase plot
	m = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:pl_rade => "Radius (Earth radii)",
	)
	m_ACCESS = mapping(
		color = :status => 
		sorter("Observing", "Data complete", "Published") => "Status",
	)
	marker_all = visual(color=(:darkgrey, 0.25))
	marker_ACCESS = visual(markersize=18)
	plt = m*(
		data(df)*marker_all +
		data(df_ACCESS) * m_ACCESS * marker_ACCESS
	)
	# colors = [
	# 	"Observing"=>COLORS[2], "Data complete"=>COLORS[1], "Published"=>COLORS[3]
	# ]
	fg = draw!(ax, plt)#, palettes=(color=colors,))
	lg = legend!(fig[1, 1], fg;
		tellwidth = false,
		halign = :right,
		valign = :center,
		framevisible = true,
		padding = (15, 15, 15, 15),
		#margin = (10, 10, 0, 0),
		#titlegap = 24,
	)
	
	# Annotate target names
	for (name, T, R) in zip(df_ACCESS.name_abbrv, df_ACCESS.pl_eqt, df_ACCESS.pl_rade) 
		if name == "W50b"
			align = (:left, :top)
		elseif name == "W96b"
			align = (:right, :baseline)
		elseif name == "W107b"
			align = (:right, :top)
		else
			align = (:left, :baseline)
		end
		text!(ax, name, position=(T, R), align=align, textsize=12)
	end
	
	fig
end

# ╔═╡ b5c0dbc8-d700-473e-9f00-d89a319f6432
md"""
## Targets with tspecs
"""

# ╔═╡ 24cd6863-9b93-473b-a66b-993e3f7d50a5
df.pl_name[1]

# ╔═╡ b4c7316d-d198-4449-ad45-66397fd1a9a5
tspec_targs = [
	"GJ 436 b",
	"GJ 1214 b",
	"GJ 3470 b",
	"HAT-P-1 b",
	"HAT-P-11 b",
	"HAT-P-12 b",
	"HAT-P-23 b",
	"HAT-P-26 b",
	"HAT-P-32 b",
	"HAT-P-38 b",
	"HAT-P-41 b",
	"HD 97658 b",
	"HD 189733 b",
	"HD 209458 b",
	"K2-18 b",
	"KELT-11 b",
	"Kepler-51 b",
	"Kepler-51 d",
	"TRAPPIST-1 b",
	"TRAPPIST-1 c",
	"TRAPPIST-1 d",
	"TRAPPIST-1 e",
	"TRAPPIST-1 f",
	"TRAPPIST-1 g",
	"WASP-6 b",
	"WASP-12 b",
	"WASP-17 b",
	"WASP-19 b",
	"WASP-21 b",
	"WASP-31 b",
	"WASP-39 b",
	"WASP-43 b",
	"WASP-50 b",
	"WASP-52 b",
	"WASP-62 b",
	"WASP-63 b",
	"WASP-67 b",
	"WASP-76 b",
	"WASP-79 b",
	"WASP-101 b",
	"WASP-107 b",
	"WASP-121 b",
	"WASP-127 b",
	"XO-1 b",
]

# ╔═╡ 0f118d0e-0eb6-4517-8378-9623337f73ca
df_tspecs = @subset df :pl_name ∈ tspec_targs

# ╔═╡ 8d519cad-8da6-409e-b486-2bc9a6008e0f
function yee!(ax, targ; al_x=:center, al_y=:bottom, df=df_HGHJs, offset=nothing)
	targg = split(targ, " (")[1]
	x, y = val.(Ref(df), Ref(targg), [:pl_eqt, :g_SI])
	text!(ax, targ, position=(x, y), align=(al_x, al_y), textsize=16, offset=offset)
end

# ╔═╡ c0f576a7-908d-4f10-86e7-cadbb7c77c09
let
	fig = Figure(resolution=(800, 700))
	ax = Axis(fig[1, 1], limits=((0, nothing), (0, 60)))
	
	p = (
		data(df) * visual(color=(:darkgrey, 0.25)) +
		data(df_tspecs) * visual(markersize=15, color=:grey)
		) *
		mapping(
			:pl_eqt => "Equilibrium temperature (K)",
			:g_SI => "Surface gravity (m/s²)",
		)

	draw!(ax, p)

	yee!(ax, "WASP-43 b (Weaver+ 2020)", offset=(0, 5))
	yee!(ax, "HAT-P-23 b (Weaver+ 2021)", al_x=:left, offset=(0, 5))
	yee!(ax, "WASP-50 b (Weaver+ in prep.)", al_x=:right, offset=(0, 5))

	hlines!(ax, 20, color=:darkgrey, linestyle=:dash)

	path = "../../ACCESS_WASP-50b/figures/pop"
	mkpath(path)
	save("$(path)/tspec_pop.png", fig)

	fig
end

# ╔═╡ Cell order:
# ╟─cd13d7f3-0ea3-4631-afd9-5f3e359000e6
# ╟─6b06701b-05e2-4284-a308-e9edeb65a648
# ╠═f396cda3-f535-4ad9-b771-7ccbd45c54f3
# ╟─4d1a7740-24c7-4cec-b788-a386bc25f836
# ╠═2702ba80-993c-4cca-bd14-0d4d6b67362b
# ╠═e0365154-d6c8-4db2-bb85-bf2536a3aa74
# ╠═9aed232f-ec74-4ec6-9ae7-06b90539833b
# ╟─4cbbb1e8-e5fb-4ab0-a7e6-7881c2dde032
# ╟─7f956b36-ce65-4e4e-afa5-9b97b9e06954
# ╟─0f9262ef-b774-45bc-bdab-46860779683d
# ╠═c98c5618-4bac-4322-b4c3-c76181f47889
# ╟─18094afc-b77f-4cae-a30c-2691d34125d8
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╠═958453c3-7993-4620-ab7f-e7ad79781dd5
# ╠═d62b5506-1411-49f2-afe3-d4aec70641a1
# ╠═f07ad06b-81d2-454f-988f-a7ae1713eac4
# ╠═2776646e-47e7-4b9e-ab91-4035bc6df99f
# ╠═c7960066-cc33-480c-807b-c56ead4262bf
# ╠═abaee9cc-9841-4b6b-ad33-2093c27422c8
# ╠═1e587a84-ed43-4fac-81cf-134a4f3d65d5
# ╠═c1e63cf3-7f30-4858-bdd6-125d2a99529f
# ╠═d2d6452b-426c-4c61-8dc4-d210c0cd9d41
# ╟─2476f26e-f7cc-4f8e-ac66-60b85a46cb2c
# ╟─8e4b1149-abf1-4ded-b20f-4765d2ee49a9
# ╠═8cece13d-34cb-40df-8986-ac0dad210e58
# ╠═df499885-3f7d-4293-bcf5-9c761ac3ccd2
# ╠═60016c3f-6968-4d3c-ac73-d324b2a071e0
# ╠═d6791747-7503-49a9-903c-c479fc0c3d49
# ╟─b5c0dbc8-d700-473e-9f00-d89a319f6432
# ╠═24cd6863-9b93-473b-a66b-993e3f7d50a5
# ╠═b4c7316d-d198-4449-ad45-66397fd1a9a5
# ╠═0f118d0e-0eb6-4517-8378-9623337f73ca
# ╠═c0f576a7-908d-4f10-86e7-cadbb7c77c09
# ╠═8d519cad-8da6-409e-b486-2bc9a6008e0f
# ╠═24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
