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
	using Unitful: k, σ, G
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
		"pl_orbper",
		"pl_eqt",
		"pl_dens",
		"pl_trandep",
		"pl_trandur",
		"sy_jmag",
		"sy_tmag",
		"st_teff",
		"st_rad",
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

# ╔═╡ 4cbbb1e8-e5fb-4ab0-a7e6-7881c2dde032
md"""
With this list of $(nrow(df_all)) transiting exoplanets with observed stellar and planetary parameters, we now turn to visualizing the subset of each population from the list above.
"""

# ╔═╡ 7c431d58-e262-4d93-af28-1f65cfe9e452
g_low, g_high = 20.0, 100.0 

# ╔═╡ 6b9906d6-eaa5-4994-a751-296c6dcf9570
TSM_low, TSM_high = 1.0, 15.0

# ╔═╡ 7f956b36-ce65-4e4e-afa5-9b97b9e06954
md"""
## HGHJ population 💪

We define this population to be the sample of transiting exoplanets with surface gravities between $(g_low) – $(g_high) m/s², with a TSM relative to HAT-P-23b between $(TSM_low) – $(TSM_high). We choose our relative TSM (TSMR) to be based on the TSM of HAT-P-23b because that is the HGHJ with the smallest TSM, that also has transmission spectra data available.
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

# ╔═╡ 46942e93-fa05-48ce-9639-248ccb63fa30
# Compute g, in m/s², assuming M, R in Earth units
compute_ge_SI(Mpe, Rpe) = Ge * Mpe / Rpe^2

# ╔═╡ 1e587a84-ed43-4fac-81cf-134a4f3d65d5
compute_H(μ, T, g) = k * T / μ #(μ * g) # For testing units

# ╔═╡ 09e280b6-7e1b-4565-9aaf-034d3420140c
compute_H(2u"u", 1_000u"K", 30u"m/s^2")# |> u"km"

# ╔═╡ 32a33cc4-e5a4-4217-8ce9-180c9d19d13f
const k_SI = ustrip(k) # Already in SI

# ╔═╡ 7056d308-670f-4d55-98a8-1be8192c1ad8
# Compute H
compute_H_SI(μ, T, g) = k_SI * T / μ #(μ * g)

# ╔═╡ f6bbfc30-21d4-4975-baf3-950af588337d
compute_H_SI(2, 1_000, 30)

# ╔═╡ 8fe94715-722b-4a3d-a9bc-4ea2fe4cc4db
u_SI = ustrip(u"kg", 1u"u")

# ╔═╡ c7960066-cc33-480c-807b-c56ead4262bf
# compute TSM, assuming M, R in Earth units, T in Kelvin
function compute_TSM(Rp, Teq, Mp, Rs, J; denom=1.0)
	f = compute_scale_factor(Rp)
	return f * (Rp^3 * Teq / (Mp * Rs^2)) * 10.0^(-J/5.0) / denom
end

# ╔═╡ 9aed232f-ec74-4ec6-9ae7-06b90539833b
begin
	df = @chain df_all begin
		dropmissing([:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag, :pl_trandep])
		@transform begin
			:pl_g = compute_g.(:pl_bmasse, :pl_rade)
			:pl_TSM = compute_TSM.(
				:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag
			)
		end
	end
	TSM0 = df.pl_TSM[df.pl_name .== "HAT-P-23 b"][1]
	df.pl_TSMR = df.pl_TSM / TSM0
	df
end

# ╔═╡ c98c5618-4bac-4322-b4c3-c76181f47889
df_HGHJs_all = @chain df begin
	@subset @. (TSM_low ≤ :pl_TSMR ≤ TSM_high) & (g_low ≤ :pl_g ≤ g_high)
	sort(:pl_TSMR, rev=true)
end

# ╔═╡ d62b5506-1411-49f2-afe3-d4aec70641a1
df_HGHJs = @subset(
	df_HGHJs_all, :pl_name .∈ Ref(["HAT-P-23 b", "WASP-43 b", "WASP-50 b"])
)

# ╔═╡ c1cd9292-28b9-4206-b128-608aaf30ff9c
let
	m = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:pl_g => "Surface gravity (m/s²)",
		color = :sy_tmag => "Tess magnitude",
		markersize = :pl_TSMR => (x -> 15.0*x),
	)
	marker_open = visual(marker='○')
	marker_closed = visual(marker='●')
	plt = m * (data(df_HGHJs_all)*marker_open + data(df_HGHJs)*marker_closed)
	fg = draw(plt)
	
	ax = fg.grid[1, 1].axis
	hlines!(ax, 30.0, color=:darkgrey, linestyle=:dash)
	TSMR_range = round.(extrema(df.pl_TSMR))
	text!(ax, "TSMR: $(TSMR_range[1]) – $(TSMR_range[2])", position=(450, 45))
	
	fg
end

# ╔═╡ 2a59198c-95aa-4360-be05-49b38f1c9171
const Ge_SI = ustrip(u"m*Rearth^2/Mearth/s^2", G)

# ╔═╡ dcb16743-4ea5-45e8-936e-c6445bfa010f
begin
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

# ╔═╡ 4bd35a25-5e17-4b71-9f42-5bb5e6caa53c
let
	#df = DataFrameMacros.titanic()
	
	df0 = DataFrame(
		:name => ["apple", "banana", "watermelon"],
		:x =>[1, 2, 3],
		:y => [2, 4, 6],
	)
	
	df = @transform df0 begin
		:z = :x .* :y
		:w = :x
	end
end

# ╔═╡ Cell order:
# ╟─cd13d7f3-0ea3-4631-afd9-5f3e359000e6
# ╟─6b06701b-05e2-4284-a308-e9edeb65a648
# ╠═f396cda3-f535-4ad9-b771-7ccbd45c54f3
# ╟─4d1a7740-24c7-4cec-b788-a386bc25f836
# ╠═9aed232f-ec74-4ec6-9ae7-06b90539833b
# ╟─4cbbb1e8-e5fb-4ab0-a7e6-7881c2dde032
# ╟─7f956b36-ce65-4e4e-afa5-9b97b9e06954
# ╠═7c431d58-e262-4d93-af28-1f65cfe9e452
# ╠═6b9906d6-eaa5-4994-a751-296c6dcf9570
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╠═c98c5618-4bac-4322-b4c3-c76181f47889
# ╠═d62b5506-1411-49f2-afe3-d4aec70641a1
# ╠═2776646e-47e7-4b9e-ab91-4035bc6df99f
# ╠═c7960066-cc33-480c-807b-c56ead4262bf
# ╠═46942e93-fa05-48ce-9639-248ccb63fa30
# ╠═7056d308-670f-4d55-98a8-1be8192c1ad8
# ╠═1e587a84-ed43-4fac-81cf-134a4f3d65d5
# ╠═09e280b6-7e1b-4565-9aaf-034d3420140c
# ╠═f6bbfc30-21d4-4975-baf3-950af588337d
# ╠═32a33cc4-e5a4-4217-8ce9-180c9d19d13f
# ╠═8fe94715-722b-4a3d-a9bc-4ea2fe4cc4db
# ╠═2a59198c-95aa-4360-be05-49b38f1c9171
# ╠═dcb16743-4ea5-45e8-936e-c6445bfa010f
# ╠═24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
# ╠═4bd35a25-5e17-4b71-9f42-5bb5e6caa53c
