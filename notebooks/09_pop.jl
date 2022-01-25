### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics, CairoMakie, CSV, DataFrames, Unitful, UnitfulAstro
	using DataFrameMacros, Chain
	using HTTP
	using Unitful
	using PlutoUI
	import MarkdownLiteral: @mdx
	using Unitful: k, G
	using NaturalSort
	using Latexify
end

# ╔═╡ 7493bb13-ee41-4798-99f6-dc1df97bd624
begin
	const DATA_DIR = "data/pop"
	const FIG_DIR = "figures/pop"
	TableOfContents()
end

# ╔═╡ cd13d7f3-0ea3-4631-afd9-5f3e359000e6
@mdx """
# HGJH population

In this notebook we will explore the possible targets amenable to atmopsheric characterization near the high-gravity hot-Jupiter (HGHJ) parameter space (Tₚ ∼ 1000 K, g ∼ 30 m/s²).

!!! note "Data download"
	```
	rclone sync -P drive_ACCESS:papers/WASP-50b/$(DATA_DIR) $(DATA_DIR)
	```
	* [Direct link](https://app.box.com/s/p8crolyu1avcbpfv49n0iehcnp8ym72p)
"""

# ╔═╡ 6b06701b-05e2-4284-a308-e9edeb65a648
@mdx """
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
@mdx """
We next compute some relevant quantities from this table to help organize each population:
"""

# ╔═╡ c7eabcc6-5139-448d-abdb-ec752788bd59
strip_u(u) = x -> ustrip(u, x)

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
@mdx """
With this list of $(nrow(df_all)) transiting exoplanets with observed stellar and planetary parameters, we now turn to visualizing the subset of each population from the list above.
"""

# ╔═╡ b5c0dbc8-d700-473e-9f00-d89a319f6432
@mdx """
## Targets with tspecs
"""

# ╔═╡ b4c7316d-d198-4449-ad45-66397fd1a9a5
tspec_targs = [
	# "GJ 436 b",
	# "GJ 1214 b",
	# "GJ 3470 b",
	# "HAT-P-1 b",
	# "HAT-P-11 b",
	# "HAT-P-12 b",
	"HAT-P-23 b",
	# "HAT-P-26 b",
	# "HAT-P-32 b",
	# "HAT-P-38 b",
	# "HAT-P-41 b",
	# "HD 97658 b",
	# "HD 189733 b",
	# "HD 209458 b",
	# "K2-18 b",
	# "KELT-11 b",
	# "Kepler-51 b",
	# "Kepler-51 d",
	# "TRAPPIST-1 b",
	# "TRAPPIST-1 c",
	# "TRAPPIST-1 d",
	# "TRAPPIST-1 e",
	# "TRAPPIST-1 f",
	# "TRAPPIST-1 g",
	# "WASP-6 b",
	# "WASP-12 b",
	# "WASP-17 b",
	# "WASP-19 b",
	# "WASP-21 b",
	# "WASP-31 b",
	# "WASP-39 b",
	#"WASP-43 b",
	"WASP-50 b",
	# "WASP-52 b",
	# "WASP-62 b",
	# "WASP-63 b",
	# "WASP-67 b",
	# "WASP-76 b",
	# "WASP-79 b",
	# "WASP-101 b",
	# "WASP-107 b",
	# "WASP-121 b",
	# "WASP-127 b",
	# "XO-1 b",
]

# ╔═╡ 5721eb48-a63a-4686-b0fa-b0b76e78ad85
 1357.98 - 1343.32

# ╔═╡ ddd8abbb-f057-4b60-bc1b-ee7f51aaa70a
df_wakeford = let
	df = CSV.read("data/pop/H2O_J_data.txt", DataFrame;
		delim = " ",
		ignorerepeated = true,
	)

	df.g_SI = @. 10^df.logg / 100

	rename!(df, :T => :pl_eqt)

	df
end

# ╔═╡ 15980829-67b5-4ccc-b914-dc188df67563
sort(df_wakeford, :Name, lt=natural)

# ╔═╡ 7f956b36-ce65-4e4e-afa5-9b97b9e06954
@mdx """
## HGHJ population

We define this population to be the sample of transiting exoplanets with ... We defined our relative TSM (TSMR) to be based on the TSM of HAT-P-23b for comparison.
"""

# ╔═╡ 0f9262ef-b774-45bc-bdab-46860779683d
@mdx """
!!! note
	Inspired from [warm-worlds](https://github.com/nespinoza/warm-worlds)
"""

# ╔═╡ 18094afc-b77f-4cae-a30c-2691d34125d8
@mdx """
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

# ╔═╡ 7336f748-5a5a-476e-80d0-cb6200aefeff
_df = @chain df_all begin
	dropmissing([:pl_rade, :pl_bmasse, :pl_orbsmax, :st_teff, :sy_jmag])
	@transform begin
		:pl_eqt = compute_Teq(
			:st_teff*u"K", :st_rad*u"Rsun", :pl_orbsmax*u"AU"; α=0.1
		) |> strip_u(u"K")

		:g_SI = compute_g(
			:pl_bmasse*u"Mearth", :pl_rade*u"Rearth"
		) |> strip_u(u"m/s^2")
	end
	@transform begin
		:H_km = compute_H(
			2.0*u"u", :pl_eqt*u"K", :g_SI*u"m/s^2"
		) |> strip_u(u"km")

		:TSM = compute_TSM(
			:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag
		)
	end
end;

# ╔═╡ c1e63cf3-7f30-4858-bdd6-125d2a99529f
compute_ΔD(H, Rₚ, Rₛ) = 2.0 * H * Rₚ/Rₛ^2

# ╔═╡ 56ad4c87-069b-4815-955b-7a8d7d012031
df = @chain _df begin
	@aside begin
		targ_idx = _.pl_name .== "HAT-P-23 b"
		TSM_target = _.TSM[targ_idx][1]
	end
	@transform begin
		:TSMR = :TSM / TSM_target

		:ΔD_ppm = compute_ΔD(
			:H_km*u"km", :pl_rade*u"Rearth", :st_rad*u"Rsun"
		) |> x -> uconvert(NoUnits, x) * 5.0 * 1e6
	end
	@aside ΔD_ppm_targ = _.ΔD_ppm[targ_idx][1]
	@transform begin
		:ΔDR_ppm = :ΔD_ppm / ΔD_ppm_targ
	end
end

# ╔═╡ c98c5618-4bac-4322-b4c3-c76181f47889
df_HGHJs_all = @chain df begin
	@subset (1.0 ≤ :TSMR) &
	(15.0 ≤ :g_SI ≤ 100.0)
	sort(:TSMR, rev=true) # To stack smaller circles on top in Figure
end

# ╔═╡ d62b5506-1411-49f2-afe3-d4aec70641a1
df_HGHJs = @subset(
	df_HGHJs_all, :pl_name .∈ Ref(
		["HAT-P-23 b", "WASP-43 b", "WASP-50 b", "HD 189733 b"]
	)
)

# ╔═╡ 8d519cad-8da6-409e-b486-2bc9a6008e0f
function yee!(ax, targ; al_x=:center, al_y=:bottom, df=df_HGHJs, offset=nothing)
	targg = split(targ, " (")[1]
	x, y = val.(Ref(df), Ref(targg), [:pl_eqt, :g_SI])
	text!(ax, targ, position=(x, y), align=(al_x, al_y), textsize=16, offset=offset)
end

# ╔═╡ 53c4fd02-7a48-47a1-9341-ede3e8d497f7
y = @chain df_HGHJs_all begin
	@select :pl_name :disc_facility :pl_eqt :g_SI :ΔD_ppm :TSMR
	@sort :ΔD_ppm
	#first(20)
end

# ╔═╡ 94a5f868-d043-4c1f-831c-17ebabd3df6c
@with_terminal begin
	println(latextabular(y, latex=false, fmt="%.2f"))
end

# ╔═╡ 0f118d0e-0eb6-4517-8378-9623337f73ca
df_tspecs = @subset df :pl_name ∈ tspec_targs

# ╔═╡ b20216ab-edf0-4f70-a9ab-3b9148e4392c
sort(df_tspecs, :pl_name, lt=natural)

# ╔═╡ 72fc2033-74cf-4e62-9c26-71016caacbea
@subset df :pl_name ∈ ["HAT-P-23 b", "WASP-43 b", "WASP-50 b"]

# ╔═╡ 683a8d85-b9a8-4eab-8a4b-e2b57d0783c0
@mdx """
## Notebook setup
"""

# ╔═╡ 95bb5b9e-0c50-48fa-bf4c-d0819c327bcc
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 81e14e30-2880-40f8-b4fa-a48a2dc34db7
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (800, 600)
	const FIG_LARGE = (1_200, 1_000)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#a6cee3",  # Cyan
			"#fdbf6f",  # Yellow
			"#ff7f00",  # Orange
			"#1f78b4",  # Blue
		]
	)
	
	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (
				xlabelsize = 18,
				ylabelsize = 18,
				topspinevisible = true,
				rightspinevisible = true,
				topspinecolor = :darkgrey,
				rightspinecolor = :darkgrey
			),
			Label = (
				textsize = 18,
				padding = (0, 10, 0, 0),
				font = AlgebraOfGraphics.firasans("Medium")
			),
			Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Medium"),),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			fontsize = 18,
			rowgap = 5,
			colgap = 5,
		)
	)
	
	COLORS

end

# ╔═╡ c0f576a7-908d-4f10-86e7-cadbb7c77c09
let
	fig = Figure(resolution=FIG_WIDE)
	ax = Axis(fig[1, 1], limits=((0, 4100), (-1, 60)))
	
	p = data(df) *
			mapping(
				:pl_eqt => "Equilibrium temperature (K)",
				:g_SI => "Surface gravity (m/s²)",
			) *
			visual(color=(:darkgrey, 0.25)) +
		data(df_wakeford) *
				mapping(
				:pl_eqt => "Equilibrium temperature (K)",
				:g_SI => "Surface gravity (m/s²)",
				color = :H2OJ => "H₂O - J",
			) *
			visual(marker=:rect, markersize=20, strokewidth=1, colormap=:cividis) +
		data(df_tspecs) *
			mapping(
				:pl_eqt => "Equilibrium temperature (K)",
				:g_SI => "Surface gravity (m/s²)",
			) *
			visual(marker='□', markersize=20)

	fg = draw!(ax, p)
	colorbar!(fig[1, 2], fg)

	yee!(ax, "WASP-43 b (Weaver+ 2020)", al_x=:left, offset=(10, 0))
	yee!(ax, "HAT-P-23 b (Weaver+ 2021)", al_x=:left, offset=(0, 8))
	yee!(ax, "WASP-50 b (this work)", al_x=:right, offset=(0, 8))
	yee!(ax, "HD 189733 b (Sing+ 2016.)", al_x=:left, offset=(0, 8))

	hlines!(ax, 20, color=:darkgrey, linestyle=:dash)

    savefig(fig, "$(FIG_DIR)/t_vs_g.png")

	fig
end

# ╔═╡ c1cd9292-28b9-4206-b128-608aaf30ff9c
# TODO: Place latitude constraints
let
	fig = Figure(resolution=FIG_WIDE)
	ax = Axis(fig[1, 1], limits=((0, 4100), (-1, 60)))

	# HGHJ g boundary
	hlines!(ax, 20.0, color=:darkgrey, linestyle=:dash)
	
	# Phase plot
	markersize_factor = 8.0
	m = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:g_SI => "Surface gravity (m/s²)",
		color = :ΔD_ppm => "ΔD (ppm)",
		markersize = :TSMR => (x -> markersize_factor*x),
	)
	plt = m*data(df_HGHJs_all) * visual(colormap=:viridis)
	fg = draw!(ax, plt)
	colorbar!(fig[1, 2], fg)
	
	# HGHJ g boundary
	# ax = fg.grid[1, 1].axis
	#hlines!(ax, 20.0, color=:darkgrey, linestyle=:dash)
	
	# Annotate HGHJs with tspec observations
	HP23x, HP23y = val.(Ref(df_HGHJs), Ref("HAT-P-23 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"HAT-P-23b",
		Point2f((HP23x[1], HP23y[1])) .- (-600, -3),
		Point2f((HP23x[1], HP23y[1])),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	W43x, W43y = val.(Ref(df_HGHJs), Ref("WASP-43 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"WASP-43b",
		Point2f((W43x[1], W43y[1])) .- (800, 5),
		Point2f((W43x[1], W43y[1])),
		0.3,
		0.1;
		align = (:center, :top),
	)
	W50x, W50y = val.(Ref(df_HGHJs), Ref("WASP-50 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"WASP-50b",
		Point2f((W50x[1], W50y[1])) .- (-300, -8),
		Point2f((W50x[1], W50y[1])),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	hdx, hdy = val.(Ref(df_HGHJs), Ref("HD 189733 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"HD 189733 b",
		Point2f((hdx[1], hdy[1])) .- (600, -4),
		Point2f((hdx[1], hdy[1])),
		0.0,
		0.0;
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
		patchsize = (120, 80),
		framevisible = true,
		padding = (5, 5, -24, 10),
		margin = (0, 0, 0, 0),
		titlegap = 24,
	)

	savefig(fig, "$(FIG_DIR)/hg_pop.png")
	
	fig
end

# ╔═╡ 1c4e56b7-0b0e-4c32-89ac-38a880bf56b9
html"""
<style>
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
# ╟─cd13d7f3-0ea3-4631-afd9-5f3e359000e6
# ╠═7493bb13-ee41-4798-99f6-dc1df97bd624
# ╟─6b06701b-05e2-4284-a308-e9edeb65a648
# ╠═f396cda3-f535-4ad9-b771-7ccbd45c54f3
# ╟─4d1a7740-24c7-4cec-b788-a386bc25f836
# ╠═7336f748-5a5a-476e-80d0-cb6200aefeff
# ╠═56ad4c87-069b-4815-955b-7a8d7d012031
# ╠═c7eabcc6-5139-448d-abdb-ec752788bd59
# ╠═c98c5618-4bac-4322-b4c3-c76181f47889
# ╠═d62b5506-1411-49f2-afe3-d4aec70641a1
# ╠═e0365154-d6c8-4db2-bb85-bf2536a3aa74
# ╠═9aed232f-ec74-4ec6-9ae7-06b90539833b
# ╟─4cbbb1e8-e5fb-4ab0-a7e6-7881c2dde032
# ╟─b5c0dbc8-d700-473e-9f00-d89a319f6432
# ╠═b4c7316d-d198-4449-ad45-66397fd1a9a5
# ╠═0f118d0e-0eb6-4517-8378-9623337f73ca
# ╠═72fc2033-74cf-4e62-9c26-71016caacbea
# ╠═5721eb48-a63a-4686-b0fa-b0b76e78ad85
# ╠═b20216ab-edf0-4f70-a9ab-3b9148e4392c
# ╠═15980829-67b5-4ccc-b914-dc188df67563
# ╠═ddd8abbb-f057-4b60-bc1b-ee7f51aaa70a
# ╠═c0f576a7-908d-4f10-86e7-cadbb7c77c09
# ╠═8d519cad-8da6-409e-b486-2bc9a6008e0f
# ╠═7f956b36-ce65-4e4e-afa5-9b97b9e06954
# ╟─0f9262ef-b774-45bc-bdab-46860779683d
# ╠═94a5f868-d043-4c1f-831c-17ebabd3df6c
# ╠═53c4fd02-7a48-47a1-9341-ede3e8d497f7
# ╟─18094afc-b77f-4cae-a30c-2691d34125d8
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╠═958453c3-7993-4620-ab7f-e7ad79781dd5
# ╠═f07ad06b-81d2-454f-988f-a7ae1713eac4
# ╠═2776646e-47e7-4b9e-ab91-4035bc6df99f
# ╠═c7960066-cc33-480c-807b-c56ead4262bf
# ╠═abaee9cc-9841-4b6b-ad33-2093c27422c8
# ╠═1e587a84-ed43-4fac-81cf-134a4f3d65d5
# ╠═c1e63cf3-7f30-4858-bdd6-125d2a99529f
# ╟─683a8d85-b9a8-4eab-8a4b-e2b57d0783c0
# ╟─95bb5b9e-0c50-48fa-bf4c-d0819c327bcc
# ╠═81e14e30-2880-40f8-b4fa-a48a2dc34db7
# ╠═24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
# ╟─1c4e56b7-0b0e-4c32-89ac-38a880bf56b9
