### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics, CairoMakie, CSV, DataFrames, Unitful, UnitfulAstro
	using DataFramesMeta
	using HTTP
	using Unitful
	using PlutoUI
	using Unitful: k, G
	
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
## Data sample 🔭

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
		"pl_eqt",
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
		df_all, [:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag]
	)
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
	#@subset @. (TSM_low ≤ :TSMR ≤ TSM_high) & (g_low ≤ :g_SI ≤ g_high)
	@subset @. (:TSMR ≥ 1.0) & (10.0 ≤ :g_SI ≤ 53) & (100 ≤ :ΔD_ppm)
	sort(:TSMR, rev=true)
end

# ╔═╡ e6d3e2b6-8895-46cd-8836-611d5cc4f5d3
extrema(df_HGHJs_all.TSMR)

# ╔═╡ 309a87d4-8d9d-4228-a2fd-f71553e129b9
describe(df_HGHJs_all)

# ╔═╡ d62b5506-1411-49f2-afe3-d4aec70641a1
df_HGHJs = @subset(
	df_HGHJs_all, :pl_name .∈ Ref(["HAT-P-23 b", "WASP-43 b", "WASP-50 b"])
)

# ╔═╡ c1cd9292-28b9-4206-b128-608aaf30ff9c
# TODO: Place latitude constraints
let
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
	plt = m*data(df_HGHJs_all) + data(df_HGHJs)*marker_closed*m2
	fg = draw(plt)
	
	# HGHJ g boundary
	ax = fg.grid[1, 1].axis
	hlines!(ax, 20.0, color=:darkgrey, linestyle=:dash)
	
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
		#margin = (10, 10, 0, 0),
		titlegap = 24,
	)
	
	fg
end

# ╔═╡ 2476f26e-f7cc-4f8e-ac66-60b85a46cb2c
md"""
## CO – CH₄ transition
"""

# ╔═╡ f7c0732f-9ae1-45f0-a51b-271197f3e576
df_CO_CH₄_all = @chain df begin
	@subset @. (300 ≤ :ΔD_ppm ≤ 1700)
	sort(:TSMR, rev=true)
end

# ╔═╡ 60016c3f-6968-4d3c-ac73-d324b2a071e0
let
	# Phase plot
	markersize_factor = 15.0
	m = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:pl_rade => "Radius (Rₑₐᵣₜₕ)",
		color = :ΔD_ppm => "ΔD (ppm)",
		markersize = :TSMR => (x -> markersize_factor*x),
	)
	marker_open, marker_closed = visual(marker='○'), visual(marker='●')
	plt = m * (data(df_CO_CH₄_all)*marker_open)# + data(df_HGHJs)*marker_closed)
	fg = draw(plt)
	
	# HGHJ g boundary
	ax = fg.grid[1, 1].axis
	vlines!(ax, 1_400.0, color=:darkgrey, linestyle=:dash)
	
	# TSMR legend
	tsmrs = [4, 2, 1]
	axislegend(
		ax,
		[MarkerElement(marker='○', markersize=markersize_factor*ms) for ms ∈ tsmrs],
		["$tsmr" for tsmr ∈ tsmrs],
		"TSMR",
		position = :lt,
		patchsize = (60, 55),
		framevisible = true,
		padding = (5, 5, 8, 24),
		margin = (10, 0, 0, 0)
	)
	
	fg
end

# ╔═╡ Cell order:
# ╟─cd13d7f3-0ea3-4631-afd9-5f3e359000e6
# ╟─6b06701b-05e2-4284-a308-e9edeb65a648
# ╠═f396cda3-f535-4ad9-b771-7ccbd45c54f3
# ╟─4d1a7740-24c7-4cec-b788-a386bc25f836
# ╠═2702ba80-993c-4cca-bd14-0d4d6b67362b
# ╠═9aed232f-ec74-4ec6-9ae7-06b90539833b
# ╟─4cbbb1e8-e5fb-4ab0-a7e6-7881c2dde032
# ╟─7f956b36-ce65-4e4e-afa5-9b97b9e06954
# ╠═e6d3e2b6-8895-46cd-8836-611d5cc4f5d3
# ╟─0f9262ef-b774-45bc-bdab-46860779683d
# ╠═309a87d4-8d9d-4228-a2fd-f71553e129b9
# ╠═c98c5618-4bac-4322-b4c3-c76181f47889
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╠═d62b5506-1411-49f2-afe3-d4aec70641a1
# ╠═2776646e-47e7-4b9e-ab91-4035bc6df99f
# ╠═c7960066-cc33-480c-807b-c56ead4262bf
# ╠═abaee9cc-9841-4b6b-ad33-2093c27422c8
# ╠═1e587a84-ed43-4fac-81cf-134a4f3d65d5
# ╠═c1e63cf3-7f30-4858-bdd6-125d2a99529f
# ╟─2476f26e-f7cc-4f8e-ac66-60b85a46cb2c
# ╠═f7c0732f-9ae1-45f0-a51b-271197f3e576
# ╠═60016c3f-6968-4d3c-ac73-d324b2a071e0
# ╠═24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
