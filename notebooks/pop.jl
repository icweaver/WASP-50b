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

# ╔═╡ c64b8c73-786c-451b-b8fc-43ee5ded5fcd
md"""
Column name definitions [here](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html)
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

# ╔═╡ 9aed232f-ec74-4ec6-9ae7-06b90539833b
df = @chain df_all begin
	dropmissing([:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag, :pl_trandep])
	#@subset(:pl_rade .≤ 11.0)
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
function compute_TSM(Rp, Teq, Mp, Rs, J; denom=1.0)
	f = compute_scale_factor(Rp)
	#f = 1.0
	return f * (Rp^3 * Teq / (Mp * Rs^2)) * 10.0^(-J/5.0) / denom
end

# ╔═╡ bb90ce6d-a354-425e-abfd-5f7e69997f22
TSM₀ = @subset(df, :pl_name .== "HAT-P-23 b").pl_TSM[1]

# ╔═╡ f37b0b30-4d62-4a72-a033-4c5b68b19e1a
df.pl_TSMR = df.pl_TSM ./ TSM₀

# ╔═╡ de1cc207-ccc9-41ba-bd73-e1f912b65a04
df_candidates = @chain df begin
	@subset @. (1.0 ≤ :pl_TSMR ≤ 20.0) & (20.0 ≤ :pl_g ≤ 100.0)
	sort(:pl_TSMR, rev=true)
end

# ╔═╡ d62b5506-1411-49f2-afe3-d4aec70641a1
df_HGHJs = @subset(df_candidates, :pl_name .∈ Ref(["HAT-P-23 b", "WASP-43 b", "WASP-50 b"]))

# ╔═╡ c1cd9292-28b9-4206-b128-608aaf30ff9c
let
	m = mapping(
		:pl_g => "Surface gravity (m/s²)",
		:pl_eqt => "Equilibrium temperature (K)",
		color = :sy_tmag => "Tess magnitude",
		markersize = :pl_TSMR => (x -> 20.0*x),
	)
	marker_open = visual(marker='○')
	marker_closed = visual(marker='●')
	plt = m * (data(df_candidates)*marker_open + data(df_HGHJs)*marker_closed)
	fg = draw(plt)
	
	ax = fg.grid[1, 1].axis
	vlines!(ax, 30.0, color=:darkgrey, linestyle=:dash)
	TSMR_range = round.(extrema(df_candidates.pl_TSMR))
	text!(ax, "TSMR: $(TSMR_range[1]) -- $(TSMR_range[2])", position=(35, 2_500))
	
	fg
end

# ╔═╡ 2a59198c-95aa-4360-be05-49b38f1c9171
const G = ustrip(u"m*Rearth^2/Mearth/s^2", Unitful.G)

# ╔═╡ 46942e93-fa05-48ce-9639-248ccb63fa30
compute_g(Mp, Rp) = G * Mp / Rp^2

# ╔═╡ 7b52e0aa-cf3e-472e-9ca6-06db018ac86d
@transform! df begin
	:pl_g = compute_g.(:pl_bmasse, :pl_rade)
	:pl_TSM = compute_TSM.(
		:pl_rade,
		:pl_eqt,
		:pl_bmasse,
		:st_rad,
		:sy_jmag,
	)
end

# ╔═╡ Cell order:
# ╟─c64b8c73-786c-451b-b8fc-43ee5ded5fcd
# ╠═f396cda3-f535-4ad9-b771-7ccbd45c54f3
# ╠═9aed232f-ec74-4ec6-9ae7-06b90539833b
# ╠═2776646e-47e7-4b9e-ab91-4035bc6df99f
# ╠═46942e93-fa05-48ce-9639-248ccb63fa30
# ╠═7b52e0aa-cf3e-472e-9ca6-06db018ac86d
# ╠═c7960066-cc33-480c-807b-c56ead4262bf
# ╠═bb90ce6d-a354-425e-abfd-5f7e69997f22
# ╠═f37b0b30-4d62-4a72-a033-4c5b68b19e1a
# ╠═de1cc207-ccc9-41ba-bd73-e1f912b65a04
# ╠═d62b5506-1411-49f2-afe3-d4aec70641a1
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╠═2a59198c-95aa-4360-be05-49b38f1c9171
# ╠═24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
