### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ 24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using PlutoUI
	import MarkdownLiteral: @mdx
	using AlgebraOfGraphics, CairoMakie
	using CSV, Chain, DataFrames, DataFramesMeta
	using HTTP
	using Unitful, UnitfulAstro
	import Unitful: k, G
	using NaturalSort
	using Latexify
end

# ╔═╡ 333d13f6-e6dd-4312-9d0f-84a41aebbde7
using Measurements

# ╔═╡ 7493bb13-ee41-4798-99f6-dc1df97bd624
begin
	const DATA_DIR = "data/pop"
	const FIG_DIR = "figures/pop"
	TableOfContents()
end

# ╔═╡ cd13d7f3-0ea3-4631-afd9-5f3e359000e6
@mdx """
# HGJH population

In this notebook we will explore possible targets in the high-gravity hot-Jupiter (HGHJ) parameter space (Tₚ ∼ 1000 K, g ∼ 30 m/s²) that are amenable to follow-up atmopsheric characterization.

!!! note "Data download"
	```
	rclone sync -P drive_ACCESS:papers/WASP-50b/$(DATA_DIR) $(DATA_DIR)
	```
	* [Direct link](https://app.box.com/s/p8crolyu1avcbpfv49n0iehcnp8ym72p)
"""

# ╔═╡ 6b06701b-05e2-4284-a308-e9edeb65a648
@mdx """
## Data sample

We draw our sample from the NASA Exoplanet Archive [TAP API](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html):
"""

# ╔═╡ f396cda3-f535-4ad9-b771-7ccbd45c54f3
df_all = let
		columns = [
		"pl_name",
		"disc_facility",
		"tic_id",
		"pl_rade",
		"pl_radeerr1",
		"pl_radeerr2",
		"pl_bmasse",
		"pl_bmasseerr1",
		"pl_bmasseerr2",
		"pl_orbsmax",
		"st_teff",
		"st_rad",
		"st_raderr1",
		"st_raderr2",
		"sy_jmag",
		"pl_eqt",
		"pl_eqterr1",
		"pl_eqterr2",
		"pl_eqt_reflink",
		"pl_ratror",
		"pl_ratrorerr1",
		"pl_ratrorerr2",
	]
	url = "https://exoplanetarchive.ipac.caltech.edu/TAP"
	#cond = "tran_flag+=1+and+pl_eqt+<+1000+and+pl_rade+<+4"
	cond = "tran_flag+=1"
	query = "select+$(join(columns, ','))+from+pscomppars+where+$(cond)&format=csv"
	request = HTTP.get("$(url)/sync?query=$(query)")
	CSV.read(request.body, DataFrame)
end

# ╔═╡ 86a99042-bb9b-43e6-87ae-d76f88b10533
# df_ps = let
# 		columns = [
# 		"pl_name",
# 		"st_rad",
# 	]
# 	url = "https://exoplanetarchive.ipac.caltech.edu/TAP"
# 	#cond = "tran_flag+=1+and+pl_eqt+<+1000+and+pl_rade+<+4"
# 	cond = "pl_name+=0"
# 	query = "select+$(join(columns, ','))+from+ps+where+$(cond)&format=csv"
# 	request = HTTP.get("$(url)/sync?query=$(query)")
# 	CSV.read(request.body, DataFrame)
# end

# ╔═╡ 4d1a7740-24c7-4cec-b788-a386bc25f836
@mdx """
We next compute the relevant quantities for estimating atmospheric observability:
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

# ╔═╡ b5c0dbc8-d700-473e-9f00-d89a319f6432
@mdx """
## Targets with published transmission spectra

!!! warning "TODO"
	Update [list](https://stellarplanet.org/science/exoplanet-transmission-spectra)
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

# ╔═╡ 05d65745-6972-41fe-8308-e5c97c85692b
get_gₚ(Mₚ, RₚRₛ, Rₛ) = G * Mₚ / (RₚRₛ^2 * Rₛ^2)

# ╔═╡ 32449eeb-a772-423c-bd52-fd6b3b8e2bba
df_T_vs_g = let
	df = DataFrame(
		pl_name = ["HAT-P-23 b", "HD 189733 b", "WASP-43 b", "WASP-50 b"],
		RₚRₛ = [0.11616±0.00081, 0.1504±0.0039, 0.1588±0.0040, 0.1390±0.0006],
		Rₛ   = [1.089±0.028,     0.765±0.019,   0.6506±0.0054, 0.843±0.031  ]u"Rsun",
		Mₚ   = [2.07±0.12,       1.166±0.052,   1.998±0.079,   1.4688±0.0920]u"Mjup",
		Rₚ   = [1.224±0.037,     1.119±0.038,   1.006±0.017,   1.166±0.043  ]u"Rjup",
		Teq  = [1951±30,         1209±11,       1426.7±8.5,    1394.84±32.70],
	)
	@chain df begin
		@aside g_SI_m = @. get_gₚ(_.Mₚ, _.RₚRₛ, _.Rₛ) .|> u"m/s^2"
		@transform! begin
			:g_SI = Measurements.value.(g_SI_m) |> ustrip
			:g_SI_err = Measurements.uncertainty.(g_SI_m) |> ustrip
			:T_K = Measurements.value.(:Teq)
			:T_K_err = Measurements.uncertainty.(:Teq)
		end
	end
end

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

# ╔═╡ 920ac6ab-2fb9-4980-8563-b52c708cb65f
df_T_vs_g

# ╔═╡ 7f956b36-ce65-4e4e-afa5-9b97b9e06954
@mdx """
## HGHJ population
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
	eps = @. lm*(p2 - p1) / hyp
	text!(ax, t, position=p1, align=align)
	lines!(ax, [p1 .+ eps, p2 .- eps], color=:darkgrey, linewidth=1)
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

# ╔═╡ 6bbdf2a5-0f97-404d-aa5c-a6681e8efa61
compute_g((565.714 ± 31.782)u"Mearth", (10.424 ± 1.009)u"Rearth") |> u"m/s^2"

# ╔═╡ 1e587a84-ed43-4fac-81cf-134a4f3d65d5
compute_H(μ, T, g) = k * T / (μ * g)

# ╔═╡ 7336f748-5a5a-476e-80d0-cb6200aefeff
_df = @chain df_all begin
	dropmissing([:pl_rade, :pl_bmasse, :pl_bmasseerr1, :pl_eqt, :pl_orbsmax, :st_teff, :sy_jmag, :pl_radeerr1, :pl_radeerr2, :pl_bmasseerr2, :tic_id,
	:pl_eqterr1, :pl_eqterr2, :pl_ratror, :pl_ratrorerr1, :pl_ratrorerr2,])
	@aside begin
		Mₚ_err_max = maximum.(eachrow(abs.(_[:, [:pl_bmasseerr1, :pl_bmasseerr2]])))
		Mₚ = _.pl_bmasse .± Mₚ_err_max
		Rₚ_err_max = maximum.(eachrow(abs.(_[:, [:pl_radeerr1, :pl_radeerr2]])))
		Rₚ = _.pl_rade .± Rₚ_err_max
		g_SI_m = @. compute_g(
			Mₚ*u"Mearth", Rₚ*u"Rearth"
		) |> strip_u(u"m/s^2")
	end
	@transform begin
		# :pl_eqt = @. compute_Teq(
		# 	:st_teff*u"K", :st_rad*u"Rsun", :pl_orbsmax*u"AU"; α=0.1
		# ) |> strip_u(u"K")

		:g_SI = Measurements.value.(g_SI_m)
		:g_SI_err = Measurements.uncertainty.(g_SI_m)
	end
	@transform begin
		:H_km = @. compute_H(
			2.0*u"u", :pl_eqt*u"K", :g_SI*u"m/s^2"
		) |> strip_u(u"km")

		:TSM = @. compute_TSM(
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

# ╔═╡ 31c05377-2ed4-4ad6-910f-5dedddcbf6dc
@mdx """
With this list of $(nrow(df)) transiting exoplanets, we now select the subset of HGHJs.

!!! warning "TODO"
	Double check constraints
"""

# ╔═╡ 0f118d0e-0eb6-4517-8378-9623337f73ca
df_tspecs = @subset df :pl_name ∈ tspec_targs

# ╔═╡ 72fc2033-74cf-4e62-9c26-71016caacbea
@subset df :pl_name ∈ ["HAT-P-23 b", "WASP-43 b", "WASP-50 b"]

# ╔═╡ c98c5618-4bac-4322-b4c3-c76181f47889
df_HGHJs_all = @chain df begin
	@subset (1.0 .≤ :TSMR) .&
	(20.0 .≤ :g_SI) .&
	(1_000.0 .≤ :pl_eqt) .&
	(10.0 .≤ :pl_rade)
	sort(:TSMR, rev=true) # To stack smaller circles on top in Figure
end

# ╔═╡ d62b5506-1411-49f2-afe3-d4aec70641a1
df_HGHJs = @subset(
	df_HGHJs_all, :pl_name .∈ Ref(
		["HAT-P-23 b", "WASP-43 b", "WASP-50 b", "HD 189733 b"]
	)
)

# ╔═╡ 8d519cad-8da6-409e-b486-2bc9a6008e0f
function label_text!(ax, targ;
	al_x=:center, al_y=:bottom, df=df_HGHJs, offset=nothing
)
	targg = split(targ, " (")[1]
	x, y = val.(Ref(df), Ref(targg), [:pl_eqt, :g_SI])
	text!(ax, targ, position=(x, y), align=(al_x, al_y), offset=offset)
end

# ╔═╡ 53c4fd02-7a48-47a1-9341-ede3e8d497f7
y = @chain df_HGHJs_all begin
	@select :pl_name :disc_facility :pl_eqt :g_SI :TSMR
	#first(10)
end

# ╔═╡ 94a5f868-d043-4c1f-831c-17ebabd3df6c
latextabular(y, latex=false, fmt="%.2f") |> PlutoUI.Text

# ╔═╡ 683a8d85-b9a8-4eab-8a4b-e2b57d0783c0
@mdx """
## Notebook setup 🔧
"""

# ╔═╡ 95bb5b9e-0c50-48fa-bf4c-d0819c327bcc
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig, pt_per_unit=1)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 81e14e30-2880-40f8-b4fa-a48a2dc34db7
begin
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = 72 .* (6, 8)
	const FIG_WIDE = 72 .* (12, 6)
	const FIG_LARGE = 72 .* (12, 12)
	const COLORS_SERIES = to_colormap(:seaborn_colorblind, 9)
	const COLORS = parse.(Makie.Colors.Colorant,
		[
			"#66C2A5",  # Green
			"#FDBF6F",  # Yellow
			"#FF7F00",  # Orange
			"#1F78B4",  # Blue
		]
	)

	set_aog_theme!()
	update_theme!(
		Theme(
			Axis = (
				xlabelsize = 24,
				ylabelsize = 24,
				topspinevisible = true,
				rightspinevisible = true,
				topspinecolor = :darkgrey,
				rightspinecolor = :darkgrey,
			),
			Label = (
				textsize = 24,
				font = AlgebraOfGraphics.firasans("Medium"),
			),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Regular"), textsize=24),
			palette = (color=COLORS, patchcolor=[(c, 0.35) for c in COLORS]),
			figure_padding = (0, 1.5, 0, 0.0),
			fontsize = 24,
			rowgap = 5,
			colgap = 5,
		)
	)

	COLORS
end

# ╔═╡ c0f576a7-908d-4f10-86e7-cadbb7c77c09
let
	p = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:g_SI => "Surface gravity (m/s²)",
	) *
	(
		data(df) * visual(color=(:darkgrey, 0.25)) +
		
		data(df_wakeford)
			* mapping(color=:H2OJ => "H₂O - J")
			* visual(
				marker=:rect, markersize=20, strokewidth=1, colormap=:cividis
			) +
		
		data(df_tspecs) * visual(marker='□', markersize=20)
	) +
	data(df_T_vs_g) * mapping(:T_K, :g_SI, :T_K_err) * visual(Errorbars, direction=:x) +
	data(df_T_vs_g) * mapping(:T_K, :g_SI, :g_SI_err) * visual(Errorbars, direction=:y)
	
	fg = draw(p;
		axis = (; limits=((0, 3_400), (-1, 55)), yticks=0:10:50),
		figure = (; resolution=FIG_LARGE),
		colorbar = (; limits=(0, 2.3)),
	)
	ax = fg.grid[1].axis
	
	label_text!(ax, "WASP-43 b (Weaver+ 2020)", al_x=:left, offset=(10, 0))
	label_text!(ax, "HAT-P-23 b (Weaver+ 2021)", al_x=:left, offset=(0, 8))
	label_text!(ax, "WASP-50 b (this work)", al_x=:right, offset=(8, 8))
	label_text!(ax, "HD 189733 b (Sing+ 2016.)", al_x=:left, offset=(0, 8))
	hl = hlines!(ax, 20, color=:darkgrey, linestyle=:dash)
	translate!(hl, 0, 0, -1) # Place behind plot markers

    savefig(fg, "$(FIG_DIR)/t_vs_g.pdf")

	fg
end

# ╔═╡ c1cd9292-28b9-4206-b128-608aaf30ff9c
# TODO: Place latitude constraints
let
	# Phase plot
	markersize_factor = 14.0
	m = mapping(
		:pl_eqt => "Equilibrium temperature (K)",
		:g_SI => "Surface gravity (m/s²)",
		#color = :ΔD_ppm => "ΔD (ppm)",
		markersize = :TSMR => (x -> markersize_factor*x),
	)
	p = m * data(df_HGHJs_all) * visual(colormap=:viridis, marker='○')

	fg = draw(p;
		axis = (; limits=((0, 3_400), (-1, 55)), yticks=0:10:50),
		figure = (; resolution=FIG_LARGE),
	)
	ax = fg.grid[1].axis

	# HGHJ g boundary
	hl = hlines!(ax, 20.0, color=:darkgrey, linestyle=:dash)
	translate!(hl, 0, 0, -1) # Place behind plot markers

	# Annotate HGHJs with tspec observations
	HP23x, HP23y = val.(Ref(df_HGHJs), Ref("HAT-P-23 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"HAT-P-23b",
		(HP23x[1], HP23y[1]) .- (-100, -5),
		(HP23x[1], HP23y[1]),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	W43x, W43y = val.(Ref(df_HGHJs), Ref("WASP-43 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"WASP-43b",
		(W43x[1], W43y[1]) .- (300, 0),
		(W43x[1], W43y[1]),
		0.0,
		0.0;
		align = (:center, :top),
	)
	W50x, W50y = val.(Ref(df_HGHJs), Ref("WASP-50 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"WASP-50b",
		(W50x[1], W50y[1]) .- (400, -5),
		(W50x[1], W50y[1]),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	hdx, hdy = val.(Ref(df_HGHJs), Ref("HD 189733 b"), [:pl_eqt, :g_SI])
	annotate_text!(
		ax,
		"HD 189733 b",
		(hdx[1], hdy[1]) .- (600, -4),
		(hdx[1], hdy[1]),
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
		position = :cb,
		patchsize = (160, 200),
		framevisible = true,
		padding = (10, 10, -20, 0),
		#margin = (0, 0, 0, 0),
		titlegap = -40,
		orientation = :horizontal,
	)

	savefig(fg, "$(FIG_DIR)/hg_pop.pdf")
	
	fg
end

# ╔═╡ Cell order:
# ╟─cd13d7f3-0ea3-4631-afd9-5f3e359000e6
# ╠═7493bb13-ee41-4798-99f6-dc1df97bd624
# ╟─6b06701b-05e2-4284-a308-e9edeb65a648
# ╠═f396cda3-f535-4ad9-b771-7ccbd45c54f3
# ╠═86a99042-bb9b-43e6-87ae-d76f88b10533
# ╟─4d1a7740-24c7-4cec-b788-a386bc25f836
# ╠═7336f748-5a5a-476e-80d0-cb6200aefeff
# ╠═333d13f6-e6dd-4312-9d0f-84a41aebbde7
# ╠═56ad4c87-069b-4815-955b-7a8d7d012031
# ╟─c7eabcc6-5139-448d-abdb-ec752788bd59
# ╟─31c05377-2ed4-4ad6-910f-5dedddcbf6dc
# ╠═d62b5506-1411-49f2-afe3-d4aec70641a1
# ╠═e0365154-d6c8-4db2-bb85-bf2536a3aa74
# ╠═9aed232f-ec74-4ec6-9ae7-06b90539833b
# ╟─b5c0dbc8-d700-473e-9f00-d89a319f6432
# ╠═b4c7316d-d198-4449-ad45-66397fd1a9a5
# ╠═0f118d0e-0eb6-4517-8378-9623337f73ca
# ╠═72fc2033-74cf-4e62-9c26-71016caacbea
# ╠═6bbdf2a5-0f97-404d-aa5c-a6681e8efa61
# ╠═05d65745-6972-41fe-8308-e5c97c85692b
# ╠═32449eeb-a772-423c-bd52-fd6b3b8e2bba
# ╠═ddd8abbb-f057-4b60-bc1b-ee7f51aaa70a
# ╠═c0f576a7-908d-4f10-86e7-cadbb7c77c09
# ╠═920ac6ab-2fb9-4980-8563-b52c708cb65f
# ╠═8d519cad-8da6-409e-b486-2bc9a6008e0f
# ╟─7f956b36-ce65-4e4e-afa5-9b97b9e06954
# ╠═c98c5618-4bac-4322-b4c3-c76181f47889
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╟─0f9262ef-b774-45bc-bdab-46860779683d
# ╠═94a5f868-d043-4c1f-831c-17ebabd3df6c
# ╠═53c4fd02-7a48-47a1-9341-ede3e8d497f7
# ╟─18094afc-b77f-4cae-a30c-2691d34125d8
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
