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
	using Statistics
	using Measurements
end

# ╔═╡ c4375025-2691-40e0-b6ca-cd0472b916cd
@mdx """
As of this writing, the planetary parameters reported in the archive have not all been updated to be self-consistent with the latest TICv8 stellar measurements taken, so we will cross-reference the targets from the archive and compute them here.
"""

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

We start by drawing our sample from the list of known exoplanets on the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/docs/TAP/usingTAP.html):
"""

# ╔═╡ f396cda3-f535-4ad9-b771-7ccbd45c54f3
df_exoarchive = let
		columns = [
		"pl_name",
		"tic_id",
		"ra",
		"dec",
		"hostname",
		"pl_rade",
		"pl_eqt",
	]
	url = "https://exoplanetarchive.ipac.caltech.edu/TAP"
	#cond = "tran_flag+=1+and+pl_eqt+<+1000+and+pl_rade+<+4"
	cond = "tran_flag+=1"
	query = "select+$(join(columns, ','))+from+pscomppars+where+$(cond)&format=csv"
	request = HTTP.get("$(url)/sync?query=$(query)")
	df = CSV.read(request.body, DataFrame)
	@chain df begin
		dropmissing!
		@rtransform! :tic_id = parse(Int, split(:tic_id)[end])
		rename!(_, :tic_id => :TIC)
	end
end

# ╔═╡ 2f9012dc-880d-4d02-9e0d-faaf6cd30766
df_exoarchive_HJ = @chain df_exoarchive begin
	@subset (1_000.0 .≤ :pl_eqt) .& (10.0 .≤ :pl_rade)
	@select :pl_name :TIC :ra :dec
end

# ╔═╡ fccd3615-7726-45ac-8d5d-c4f1fc4d5425
@mdx """
!!! note
	We only use the eqt and rade as an approximate cut-off, we go back re-compute these with the updated stellar params from TICv8 later
"""

# ╔═╡ cc5d072f-c4a8-4e4e-8e61-be21931026b1
CSV.write("data/pop/exoarchive_HJ_coords.txt", df_exoarchive_HJ[:, [:ra, :dec]];
	delim = '\t',
	header = false,
)

# ╔═╡ 9de3a7bc-29c0-455c-8b9a-f5d2031838ab
df_exoachive_TICv8_HJ = let
	df = CSV.read("data/pop/exoarchive_TICv8_HJ.txt", DataFrame;
		comment = "#",
		delim = ' ',
		ignorerepeated = true,
	)
	@select! df begin
		:TIC
		:Rₛ = :Rstar .± :Rstar_err
		:Mₛ = :Mstar .± :M_star_err
		:ρₛ = :rho_star .± :rho_star_err
		:T_eff = :Teff .± :Teff_err
	end
	#@transform! df :Mₛ = @. :ρₛ * :Rₛ^3 # Already in solar units
end

# ╔═╡ 58a87832-d4fe-4dee-86e2-394c75222088
CSV.read("data/pop/exoarchive_TICv8_HJ.txt", DataFrame;
	comment = "#",
	delim = ' ',
	ignorerepeated = true,
)

# ╔═╡ 5c7a8b0a-02f0-465f-bd81-04dc1b76b875
df_HJ = leftjoin(df_exoachive_TICv8_HJ, df_exoarchive_HJ, on=:TIC)

# ╔═╡ 86a99042-bb9b-43e6-87ae-d76f88b10533
df_ps_all = let
	columns = [
		# TIC ID
		"tic_id",
		
		# P
		"pl_orbper",
		"pl_orbpererr1",
		"pl_orbpererr2",
	
		# K
		"pl_rvamp",
		"pl_rvamperr1",
		"pl_rvamperr2",
	
		# i
		"pl_orbincl",
		"pl_orbinclerr1",
		"pl_orbinclerr2",
		
		# Rₚ/Rₛ
		"pl_ratror",
		"pl_ratrorerr1",
		"pl_ratrorerr2",
	
		# References
		"pl_refname",
		"st_refname",
		"pl_pubdate",
	]
	url = "https://exoplanetarchive.ipac.caltech.edu/TAP"
	#cond = "tran_flag+=1+and+pl_eqt+<+1000+and+pl_rade+<+4"
	#cond = "tran_flag+=1+and+st_refname+LIKE+'%TICv8%'"
	cond = "tran_flag+=1"
	query = "select+$(join(columns, ','))+from+ps+where+$(cond)&format=csv"
	request = HTTP.get("$(url)/sync?query=$(query)")
	
	df = CSV.read(request.body, DataFrame)
	
	@chain df begin
		dropmissing!(:tic_id)
		@rtransform! :tic_id = parse(Int, split(:tic_id)[end])
		rename!(_, :tic_id => :TIC)
	end
end

# ╔═╡ 92fbb7d7-9782-44d4-b1b7-6db89d78a032
df_ps = leftjoin(df_HJ, df_ps_all, on=:TIC)

# ╔═╡ e8a13c3b-819a-490e-a967-e2da54ca6617
# for df in groupby(df_ps_all, :pl_name)
# 	df_rv = @chain df begin
# 		@subset :pl_pubdate .== maximum(:pl_pubdate)
# 		@select :pl_name :pl_rvamp :pl_rvamperr1 :pl_rvamperr2 :pl_refname
# 	end
# 	K, K_err = max_m(df_rv.pl_rvamp[1], df_rv.pl_rvamperr1[1], df_rv.pl_rvamperr2[1])
# 	# @aside begin
# 	# 	@subset _.pl_pubdate .== maximum(_.pl_pubdate)
# 	# end
# 	#@select :pl_name #:pl_rvamp :pl_rvamperr1 :pl_rvamperr2 :pl_refname
# 	# 	end
# 	# 	K, K_err = max_m(
# 	# 		df_rv.pl_rvamp[1], df_rv.pl_rvamperr1[1], df_rv.pl_rvamperr2[1]
# 	# 	)
# 	# end
# 	#@combine :K = K
# end

# ╔═╡ d3f2eace-6272-4281-95fd-4b931c2ae332
ρ_sun = uconvert(u"g/cm^3", (1*u"Msun" / ((4/3) * π * u"Rsun"^3))) |> ustrip

# ╔═╡ 97c9f1ae-21da-4f99-94ff-a8adaabf30bb
@mdx """
## Radial velocity (RV params)
* K
"""

# ╔═╡ 3ed05d01-b489-46d6-bcd4-9028d756ab35
@mdx """
## Transit params
* period
* inclination
* transit depth
"""

# ╔═╡ d6598eab-33b3-4873-b6fe-b16c6d5c37d7
max_m(p, pu, pd) = p ± mean(skipmissing((pu, abs(pd))))

# ╔═╡ 2584860a-8e24-49f7-a7d5-4c99c8deda8e
function extract_K(df₀)
	df = @subset df₀ :pl_pubdate .== maximum(:pl_pubdate)
	return max_m(df.pl_rvamp[1], df.pl_rvamperr1[1], df.pl_rvamperr2[1])
end

# ╔═╡ 759b0ca7-ade4-4929-afa5-51e0ab133a5b
begin
	TICs = []
	K_mps = []
	for df ∈ groupby(dropmissing(df_ps, [:pl_rvamp, :pl_rvamperr1, :pl_rvamperr2, :pl_refname, :pl_pubdate]), :TIC;
	sort = true,
	)
		K = extract_K(df)
		push!(TICs, df.TIC[1])
		push!(K_mps, K)
	end
end

# ╔═╡ 2ae250e2-cc4d-4824-9795-a8bc0a4b469b
df_K = DataFrame(; TIC=TICs, K_mps)

# ╔═╡ 0aa8aaf2-5343-4b6e-a47b-cf0fc8d27643
function extract_orb_params(df₀)
	df = @subset df₀ :pl_pubdate .== maximum(:pl_pubdate)
	return (
		max_m(df.pl_ratror[1], df.pl_ratrorerr1[1], df.pl_ratrorerr2[1]),
		max_m(df.pl_orbper[1], df.pl_orbpererr1[1], df.pl_orbpererr2[1]),
		max_m(df.pl_orbincl[1], df.pl_orbinclerr1[1], df.pl_orbinclerr2[1])
	)
end

# ╔═╡ 2d63caf3-dd64-483d-8c85-5085d7aad2ac
begin
	TICS, rs, Ps, is = [], [], [], []
	for df ∈ groupby(dropmissing(df_ps, [:pl_ratror, :pl_ratrorerr1, :pl_orbper, :pl_orbpererr1, :pl_orbincl, :pl_orbinclerr1, :pl_pubdate]), :TIC;
	sort = true,
	)
		r, P, i = extract_orb_params(df)
		push!(TICS, df.TIC[1])
		push!(rs, r)
		push!(Ps, P)
		push!(is, i)
	end
end

# ╔═╡ 7898ea90-4863-4beb-995e-7a251235aa88
df_orb = DataFrame(; TIC=TICS, r=rs, P_d=Ps, i_deg=is)

# ╔═╡ 0128ada1-ca70-4c74-8d23-0704a6ed2f77
nrow.((df_HJ, df_K, df_orb))

# ╔═╡ 4d1a7740-24c7-4cec-b788-a386bc25f836
@mdx """
We next compute the relevant quantities for estimating atmospheric observability:
"""

# ╔═╡ 542c59fd-782f-4e15-ab6c-a450bf4714ba
compute_Mₚ(K, i, P, Mₛ) = (K/sin(i)) * cbrt(P / (2.0π*G)) * Mₛ^(2//3)

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

# ╔═╡ 31ecaea8-8fe4-49d1-b7ab-11988981b7c9
@with_terminal PrettyTables.pretty_table(df_T_vs_g, nosubheader=true)

# ╔═╡ 84a9837c-552c-4982-a55c-8d268e463758
df_T_vs_g[df_T_vs_g.pl_name .== "HAT-P-23 b", :].Mₚ[1] |> u"Mearth"

# ╔═╡ ddd8abbb-f057-4b60-bc1b-ee7f51aaa70a
df_wakeford = let
	df = CSV.read("data/pop/H2O_J_data.csv", DataFrame;
		stripwhitespace = true,
	)

	#df.g_SI = @. 10^df.logg / 100

	#rename!(df, :T => :pl_eqt)

	sort!(df, lt=natural)
end

# ╔═╡ d0e6c6d7-c53f-449f-b49d-de2c22971bb7
df_wakeford_ps = sort!(leftjoin(df_wakeford, df_ps_all, on=:pl_name), :pl_name;
	lt = natural
)

# ╔═╡ aab7f672-a312-44d2-b8cc-937f877a662d
10^2.78 / 100

# ╔═╡ ad07751a-8334-4490-82b7-5cf5a1375aad
10^4.33 / 100

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

# ╔═╡ cff3a9c4-2a9f-4769-b309-f215bfff27f7
compute_gₚ(Mₚ, RₚRₛ, Rₛ) = G * Mₚ / (RₚRₛ^2 * Rₛ^2)

# ╔═╡ 463341ab-318d-402f-9545-b44cf19a75ea
df_HJ_complete = let
	df = leftjoin(leftjoin(df_orb, df_K; on=:TIC), df_HJ;
		on = :TIC
	) |> dropmissing
	@transform! df begin
		:Mₚ_J = @. compute_Mₚ(
			:K_mps*u"m/s", :i_deg*u"°", :P_d*u"d", :Mₛ*u"Msun"
		) |> u"Mjup" |> ustrip
	end
	@transform! df begin
		:gₚ_SI = @. compute_gₚ(
			:Mₚ_J*u"Mjup", :r, :Rₛ*u"Rsun"
		) |> u"m/s^2" |> ustrip
	end
end

# ╔═╡ bdc13e42-b250-4e32-a0ce-228a69b855aa
names(df_HJ_complete)

# ╔═╡ fede3db9-ff7f-4966-b76b-c9953e7e3de2
let
	df = df_HJ_complete
	df[occursin.("WASP-43", df.pl_name), :]
end

# ╔═╡ 64e896bd-5213-4f93-8e4c-e506010e27dc
let
	K, i, P = df_HJ_complete[1, [:K_mps, :i_deg, :P_d]]
end

# ╔═╡ 808c3bb1-4215-46aa-87fb-80a51cc44485
compute_gₚ(2u"Msun", 0.001, 1u"Rsun") |> u"m/s^2" |> ustrip

# ╔═╡ 1e587a84-ed43-4fac-81cf-134a4f3d65d5
compute_H(μ, T, g) = k * T / (μ * g)

# ╔═╡ 7336f748-5a5a-476e-80d0-cb6200aefeff
_df = @chain df_all begin
	dropmissing([:pl_rade, :pl_bmasse, :pl_eqt, :pl_orbsmax, :st_teff, :sy_jmag, :tic_id, :pl_ratror, :st_rad,
	:pl_bmasseerr1, :pl_bmasseerr2, :st_raderr1, :pl_eqterr1, :pl_eqterr2])
	@aside begin
		Mₚ_err_max = maximum.(
			skipmissing.(eachrow(abs.(_[:, [:pl_bmasseerr1, :pl_bmasseerr2]]))))
		Mₚ = _.pl_bmasse .± Mₚ_err_max
		Rₚ_err_max = maximum.(
			skipmissing.(eachrow(abs.(_[:, [:pl_radeerr1, :pl_radeerr2]]))))
		Rₚ = _.pl_rade .± Rₚ_err_max
		Rₛ_err_max = maximum.(
			skipmissing.(eachrow(abs.(_[:, [:st_raderr1, :st_raderr2]]))))
		Rₛ = _.st_rad .± Rₛ_err_max
		RₚRₛ_err_max = maximum.(
			skipmissing.(eachrow(abs.(_[:, [:pl_ratrorerr1, :pl_ratrorerr2]]))))
		RₚRₛ = _.pl_ratror .± RₚRₛ_err_max
		g_SI_m = @. get_gₚ(Mₚ*u"Mearth", RₚRₛ, Rₛ*u"Rsun") |> strip_u(u"m/s^2")
		T_eq_err_max = maximum.(
			skipmissing.(eachrow(abs.(_[:, [:pl_eqterr1, :pl_eqterr2]]))))
		T_eq_m = _.pl_eqt .± T_eq_err_max
	end
	@transform begin
		:pl_eqt = @. compute_Teq(
			:st_teff*u"K", :st_rad*u"Rsun", :pl_orbsmax*u"AU"; α=0.1
		) |> strip_u(u"K")
		:g_SI = Measurements.value.(g_SI_m)
		:g_SI_err = Measurements.uncertainty.(g_SI_m)
		:Teq_K = Measurements.value.(T_eq_m)
		:Teq_K_err = Measurements.uncertainty.(T_eq_m)
	end
	@transform begin
		:H_km = @. compute_H(
			2.0*u"u", :pl_eqt*u"K", :g_SI*u"m/s^2"
		) |> strip_u(u"km")

		:TSM = @. compute_TSM(
			:pl_rade, :pl_eqt, :pl_bmasse, :st_rad, :sy_jmag
		)
	end
end

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

# ╔═╡ 8a43e9a3-f9a1-4f15-a914-2b7b4f4bf3cd
df[df.pl_name .== "HAT-P-23 b", :]

# ╔═╡ 3a4112a5-f211-4160-a6e8-688d04751a42
names(df)

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

# ╔═╡ 32c0b1c9-0e99-4c43-a473-7122af2e36ab
sort(df_HGHJs_all, :pl_name)#[df_HGHJs_all.pl_name .== "WASP-140 b", :]
# check for WASP-95 b

# ╔═╡ 53c4fd02-7a48-47a1-9341-ede3e8d497f7
y = @chain df_HGHJs_all begin
	@select :pl_name :pl_eqt :g_SI :TSMR
	#first(10)
end

# ╔═╡ 94a5f868-d043-4c1f-831c-17ebabd3df6c
latextabular(y, latex=false, fmt="%.2f") |> PlutoUI.Text

# ╔═╡ 18e6c65f-b999-41d6-81d2-7398bc05ab47
@with_terminal PrettyTables.pretty_table(y, nosubheader=true)

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
	) #+
	# data(df_T_vs_g) * mapping(:T_K, :g_SI, :T_K_err) * visual(Errorbars, direction=:x) +
	# data(df_T_vs_g) * mapping(:T_K, :g_SI, :g_SI_err) * visual(Errorbars, direction=:y)
	# data(df_tspecs) * mapping(:Teq_K, :g_SI, :Teq_K_err) * visual(Errorbars, direction=:x) +
	# data(df_tspecs) * mapping(:Teq_K, :g_SI, :g_SI_err) * visual(Errorbars, direction=:y)
	
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
# ╟─c4375025-2691-40e0-b6ca-cd0472b916cd
# ╠═7493bb13-ee41-4798-99f6-dc1df97bd624
# ╟─6b06701b-05e2-4284-a308-e9edeb65a648
# ╠═f396cda3-f535-4ad9-b771-7ccbd45c54f3
# ╠═2f9012dc-880d-4d02-9e0d-faaf6cd30766
# ╟─fccd3615-7726-45ac-8d5d-c4f1fc4d5425
# ╠═cc5d072f-c4a8-4e4e-8e61-be21931026b1
# ╠═9de3a7bc-29c0-455c-8b9a-f5d2031838ab
# ╠═58a87832-d4fe-4dee-86e2-394c75222088
# ╠═5c7a8b0a-02f0-465f-bd81-04dc1b76b875
# ╠═86a99042-bb9b-43e6-87ae-d76f88b10533
# ╠═92fbb7d7-9782-44d4-b1b7-6db89d78a032
# ╠═e8a13c3b-819a-490e-a967-e2da54ca6617
# ╠═d3f2eace-6272-4281-95fd-4b931c2ae332
# ╟─97c9f1ae-21da-4f99-94ff-a8adaabf30bb
# ╠═2584860a-8e24-49f7-a7d5-4c99c8deda8e
# ╠═759b0ca7-ade4-4929-afa5-51e0ab133a5b
# ╠═2ae250e2-cc4d-4824-9795-a8bc0a4b469b
# ╟─3ed05d01-b489-46d6-bcd4-9028d756ab35
# ╠═0aa8aaf2-5343-4b6e-a47b-cf0fc8d27643
# ╠═2d63caf3-dd64-483d-8c85-5085d7aad2ac
# ╠═7898ea90-4863-4beb-995e-7a251235aa88
# ╠═0128ada1-ca70-4c74-8d23-0704a6ed2f77
# ╠═463341ab-318d-402f-9545-b44cf19a75ea
# ╠═bdc13e42-b250-4e32-a0ce-228a69b855aa
# ╠═808c3bb1-4215-46aa-87fb-80a51cc44485
# ╠═fede3db9-ff7f-4966-b76b-c9953e7e3de2
# ╠═d6598eab-33b3-4873-b6fe-b16c6d5c37d7
# ╟─4d1a7740-24c7-4cec-b788-a386bc25f836
# ╠═542c59fd-782f-4e15-ab6c-a450bf4714ba
# ╠═64e896bd-5213-4f93-8e4c-e506010e27dc
# ╠═7336f748-5a5a-476e-80d0-cb6200aefeff
# ╠═56ad4c87-069b-4815-955b-7a8d7d012031
# ╟─c7eabcc6-5139-448d-abdb-ec752788bd59
# ╟─31c05377-2ed4-4ad6-910f-5dedddcbf6dc
# ╠═d62b5506-1411-49f2-afe3-d4aec70641a1
# ╠═e0365154-d6c8-4db2-bb85-bf2536a3aa74
# ╠═9aed232f-ec74-4ec6-9ae7-06b90539833b
# ╟─b5c0dbc8-d700-473e-9f00-d89a319f6432
# ╠═d0e6c6d7-c53f-449f-b49d-de2c22971bb7
# ╠═b4c7316d-d198-4449-ad45-66397fd1a9a5
# ╠═0f118d0e-0eb6-4517-8378-9623337f73ca
# ╠═05d65745-6972-41fe-8308-e5c97c85692b
# ╠═32449eeb-a772-423c-bd52-fd6b3b8e2bba
# ╠═31ecaea8-8fe4-49d1-b7ab-11988981b7c9
# ╠═84a9837c-552c-4982-a55c-8d268e463758
# ╠═8a43e9a3-f9a1-4f15-a914-2b7b4f4bf3cd
# ╠═ddd8abbb-f057-4b60-bc1b-ee7f51aaa70a
# ╠═aab7f672-a312-44d2-b8cc-937f877a662d
# ╠═ad07751a-8334-4490-82b7-5cf5a1375aad
# ╠═c0f576a7-908d-4f10-86e7-cadbb7c77c09
# ╠═3a4112a5-f211-4160-a6e8-688d04751a42
# ╠═8d519cad-8da6-409e-b486-2bc9a6008e0f
# ╟─7f956b36-ce65-4e4e-afa5-9b97b9e06954
# ╠═c98c5618-4bac-4322-b4c3-c76181f47889
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╟─0f9262ef-b774-45bc-bdab-46860779683d
# ╠═94a5f868-d043-4c1f-831c-17ebabd3df6c
# ╠═32c0b1c9-0e99-4c43-a473-7122af2e36ab
# ╠═53c4fd02-7a48-47a1-9341-ede3e8d497f7
# ╠═18e6c65f-b999-41d6-81d2-7398bc05ab47
# ╟─18094afc-b77f-4cae-a30c-2691d34125d8
# ╠═958453c3-7993-4620-ab7f-e7ad79781dd5
# ╠═f07ad06b-81d2-454f-988f-a7ae1713eac4
# ╠═2776646e-47e7-4b9e-ab91-4035bc6df99f
# ╠═c7960066-cc33-480c-807b-c56ead4262bf
# ╠═cff3a9c4-2a9f-4769-b309-f215bfff27f7
# ╠═1e587a84-ed43-4fac-81cf-134a4f3d65d5
# ╠═c1e63cf3-7f30-4858-bdd6-125d2a99529f
# ╟─683a8d85-b9a8-4eab-8a4b-e2b57d0783c0
# ╟─95bb5b9e-0c50-48fa-bf4c-d0819c327bcc
# ╠═81e14e30-2880-40f8-b4fa-a48a2dc34db7
# ╠═24c6a2d0-0aea-11ec-2cd4-3de7ec08b83e
