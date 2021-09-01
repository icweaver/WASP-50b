### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics, CairoMakie, CSV, DataFrames, Unitful, UnitfulAstro
	using DataFramesMeta
	using PhysicalConstants.CODATA2018: G, k_B, m_u, σ
	
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

# ╔═╡ a5559bbb-11d9-4976-bece-9aaebea7a68d
df = CSV.read("data/pop/PSCompPars_2021.08.31_23.50.02.csv", DataFrame;
	comment = "#",
)

# ╔═╡ f8cf4024-660f-4ad0-8ea9-a67fa8b9ca22
df_cleaned = dropmissing(
	df,
	[:pl_radj, :pl_eqt, :pl_bmassj, :st_rad, :sy_jmag]
)

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
		f = 0.0 # Only goes up to 1O R_earth in Kempton+2018
	end
	return f
end

# ╔═╡ 0df5cd17-d6c7-47ad-8e51-5d4aaae2ec61
function compute_TSM(f, Rp, Teq, Mp, Rs, J)
	f = compute_scale_factor(Rp)
	return f * (Rp^3 * Teq / (Mp * Rs^2)) * 10.0^(-J/5.0)
end

# ╔═╡ c55f348b-280a-4a3f-bba5-85d6eb160b03
const G_RpMp = G |> u"Rearth^3/Mearth/s^2" |> ustrip

# ╔═╡ 794df6de-c602-4b69-9502-8225c45cd156
const G_SI = G |> u"m*Rearth^2/Mearth/s^2" |> ustrip

# ╔═╡ ad322153-a033-4c71-ae28-ae9f3e18fa47
const k_B_earth = k_B |> u"kg*m^2/s^2/K" |> ustrip

# ╔═╡ 40b9adc9-94f9-4426-ad9b-839175d318b7
const mp = 2.0 * m_u |> u"kg" |> ustrip # X amu

# ╔═╡ 2a32d2af-3acd-4015-b937-cc55e5dcef11
# μ ≡ 2.0 amu (in Mearth), gp in m/s^2
# returns km
compute_H(Tp, gp) = 1e-3 * k_B_earth * Tp / (mp * gp)

# ╔═╡ 46942e93-fa05-48ce-9639-248ccb63fa30
compute_g(Mp, Rp) = G_SI * Mp / Rp^2

# ╔═╡ 8f9634d6-7600-4279-ae9c-23c91ed6cd81
df_cleaned.pl_g = compute_g.(df_cleaned.pl_bmasse, df_cleaned.pl_rade)

# ╔═╡ 88dc35de-1b8e-4b0f-acb3-8cbb77d47faa
df_cleaned.pl_H = compute_H.(df_cleaned.pl_eqt, df_cleaned.pl_g)

# ╔═╡ 3c0f3bbb-6455-4867-a069-e7816ba8eba2
let
	pop = data(@subset(df_cleaned, 1.0 .≤ :pl_H .≤ 200.0)) *
		mapping(
		:pl_g => "Surface gravity (m/s²)",
		:pl_eqt => "Equilibrium temperature (K)",
		color=:pl_H => "Scale height (km)",
	)
	
	fig = draw(pop; axis=(;limits=((0, 60), (0, 3_000))))
		
	scatter!(
		fig.figure[1, 1],
		# WASP-43, HAT-P-23, WASP-50
		[47.4, 29.2, 31.61], [1440.0, 2027.0, 1381.0];
		color=:red,
		markersize=40,
		marker='x',
	)
	
	fig
end

# ╔═╡ Cell order:
# ╠═a5559bbb-11d9-4976-bece-9aaebea7a68d
# ╠═f8cf4024-660f-4ad0-8ea9-a67fa8b9ca22
# ╠═0df5cd17-d6c7-47ad-8e51-5d4aaae2ec61
# ╠═2776646e-47e7-4b9e-ab91-4035bc6df99f
# ╠═c55f348b-280a-4a3f-bba5-85d6eb160b03
# ╠═794df6de-c602-4b69-9502-8225c45cd156
# ╠═ad322153-a033-4c71-ae28-ae9f3e18fa47
# ╠═40b9adc9-94f9-4426-ad9b-839175d318b7
# ╠═2a32d2af-3acd-4015-b937-cc55e5dcef11
# ╠═46942e93-fa05-48ce-9639-248ccb63fa30
# ╠═8f9634d6-7600-4279-ae9c-23c91ed6cd81
# ╠═88dc35de-1b8e-4b0f-acb3-8cbb77d47faa
# ╠═3c0f3bbb-6455-4867-a069-e7816ba8eba2
# ╠═24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
