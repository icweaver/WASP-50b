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

# ╔═╡ 4e6c0113-255f-46fb-84ad-f8f30b79a219
using Measurements: value, uncertainty

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
		"hostname",
		"pl_name",
		"tic_id",
		"ra",
		"dec",
		"sy_jmag",
		# "pl_rade",
		# "pl_eqt",
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
# df_exoarchive_HJ = @chain df_exoarchive begin
# 	@subset (1_000.0 .≤ :pl_eqt) .& (10.0 .≤ :pl_rade)
# 	@select :pl_name :TIC :ra :dec
# end

# ╔═╡ fccd3615-7726-45ac-8d5d-c4f1fc4d5425
@mdx """
!!! note
	~~~We only use the eqt and rade as an approximate cut-off~~~, we go back re-compute these with the updated stellar params from TICv8 later
"""

# ╔═╡ 0d629db3-7370-406f-989b-7a2caca020dc
CSV.write("data/pop/exoarchive_coords.txt", df_exoarchive[:, [:ra, :dec]];
	delim = '\t',
	header = false,
)

# ╔═╡ f8b4def8-a46f-4cbc-83c0-ff44a39c1571
ρ_sun = inv((4/3)*π)

# ╔═╡ 380d05a4-35e9-4db4-b35a-b03de9e695ee
df_exoarchive_TICv8 = let
	df = CSV.read("data/pop/exoarchive_TICv8.txt", DataFrame;
		comment = "#",
		delim = ' ',
		ignorerepeated = true,
	)
	
	df = leftjoin(df_exoarchive, df, on=:TIC) |> dropmissing! |> unique!
	
	@select! df begin
		:TIC
		:pl_name
		:st_rad = (:Rstar .± :Rstar_err)u"Rsun"
		:st_mass = (:Mstar .± :M_star_err)u"Msun"
		:st_dens_solar = ρ_sun .* (:rho_star .± :rho_star_err)u"Msun/Rsun^3"
		:st_teff = (:Teff .± :Teff_err)u"K"
		:sy_jmag
	end
	#@transform! df :Mₛ = @. :ρₛ * :Rₛ^3 # Already in solar units
end

# ╔═╡ 86a99042-bb9b-43e6-87ae-d76f88b10533
df_ps_all = let
	columns = [
		# TIC ID
		"tic_id",
		"pl_name",
		
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
		#Rₚ, Rₛ (manual fallback for Wakeford+ targs with missing data)
		"pl_radj",
		"st_rad",
		"pl_trandep",
	
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
df_ps = leftjoin(df_ps_all, df_exoarchive_TICv8, on=[:TIC, :pl_name];
makeunique=true)# |> dropmissing

# ╔═╡ 97c9f1ae-21da-4f99-94ff-a8adaabf30bb
@mdx """
## Radial velocity (RV params)
* K
"""

# ╔═╡ 2584860a-8e24-49f7-a7d5-4c99c8deda8e
function extract_K(df₀)
	df = @subset df₀ :pl_pubdate .== maximum(:pl_pubdate)
	return df.pl_rvamp[1]u"m/s", df.pl_refname[1], df.st_refname[1]
	#return max_m(df.pl_rvamp[1], df.pl_rvamperr1[1], df.pl_rvamperr2[1])
end

# ╔═╡ 759b0ca7-ade4-4929-afa5-51e0ab133a5b
begin
	pl_names_K, pl_refnames_K, st_refnames_K = String[], String[], String[]
	K_mps = []
	for df ∈ groupby(dropmissing(df_ps, [:pl_rvamp, :pl_refname, :pl_pubdate]), 		:pl_name;
		sort = true,
	)
	# for df ∈ groupby(dropmissing(df_ps, [:pl_rvamp, :pl_rvamperr1, :pl_rvamperr2, :pl_refname, :pl_pubdate]), :pl_name;
	# sort = true,
	# )
		K, pl_refname_K, st_refname_K = extract_K(df)
		push!(pl_names_K, df.pl_name[1])
		push!(K_mps, K)
		push!(pl_refnames_K, pl_refname_K)
		push!(st_refnames_K, st_refname_K)
	end
end

# ╔═╡ 2ae250e2-cc4d-4824-9795-a8bc0a4b469b
# Fixed params
df_K = DataFrame(; pl_name=pl_names_K, pl_rvamp=K_mps, pl_refnames_K, st_refnames_K)

# ╔═╡ 3ed05d01-b489-46d6-bcd4-9028d756ab35
@mdx """
## Transit params
* period
* inclination
* transit depth
"""

# ╔═╡ 0aa8aaf2-5343-4b6e-a47b-cf0fc8d27643
function extract_orb_params(df₀)
	df = @subset df₀ :pl_pubdate .== maximum(:pl_pubdate)
	
	return (df.pl_ratror[1], df.pl_orbper[1]u"d", df.pl_orbincl[1]u"°",
		df.pl_refname[1], df.st_refname[1])
	# return (
	# 	max_m(df.pl_ratror[1], df.pl_ratrorerr1[1], df.pl_ratrorerr2[1]),
	# 	max_m(df.pl_orbper[1], df.pl_orbpererr1[1], df.pl_orbpererr2[1]),
	# 	max_m(df.pl_orbincl[1], df.pl_orbinclerr1[1], df.pl_orbinclerr2[1])
	# )
end

# ╔═╡ 0c793036-d7c6-4a56-9ef6-f58b02e6530c
function extract_orb_params_wakeford(df₀)
	df = @subset df₀ :pl_pubdate .== maximum(:pl_pubdate)
	δ = df.pl_trandep[1]
	r = sqrt(δ)
	return r, df.pl_orbper[1]u"d", df.pl_orbincl[1]u"°", df.pl_refname[1], df.st_refname[1]
end

# ╔═╡ 2d63caf3-dd64-483d-8c85-5085d7aad2ac
begin
	pl_names_orb, rs, Ps, i_degs = String[], Float64[], [], []
	pl_refnames_orb, st_refnames_orb = String[], String[]
	
	for df ∈ groupby(dropmissing(df_ps,
		[:pl_ratror, :pl_orbper, :pl_orbincl, :pl_pubdate]), :TIC;
		sort = true,
	)
	# for df ∈ groupby(dropmissing(df_ps, [:pl_ratror, :pl_ratrorerr1, :pl_orbper, :pl_orbpererr1, :pl_orbincl, :pl_orbinclerr1, :pl_pubdate]), :TIC;
	# sort = true,
	# )
		r, P, i_deg, pl_refname_orb, st_refname_orb = extract_orb_params(df)
		push!(pl_names_orb, df.pl_name[1])
		push!(rs, r)
		push!(Ps, P)
		push!(i_degs, i_deg)
		push!(pl_refnames_orb, pl_refname_orb)
		push!(st_refnames_orb, st_refname_orb)
	end
end

# ╔═╡ f306bc5a-5597-4a52-b6e3-364a4230367d
function wakeford_orb(df_ps)
	pl_names, rs, Ps, i_degs = String[], Float64[], [], []
	pl_refnames, st_refnames = String[], String[]
	
	for df ∈ groupby(dropmissing(df_ps,
		[:pl_trandep, :pl_orbper, :pl_orbincl, :pl_pubdate]), :TIC;
		sort = true,
	)
		if df.pl_name[1] ∈ [
			"WASP-31 b", "WASP-63 b",
			"WASP-69 b", "WASP-76 b",
			"WASP-79 b", "WASP-101 b"
		]
			r, P, i_deg, pl_refname, st_refname = extract_orb_params_wakeford(df)
			push!(pl_names, df.pl_name[1])
			push!(rs, r)
			push!(Ps, P)
			push!(i_degs, i_deg)
			push!(pl_refnames, pl_refname)
			push!(st_refnames, st_refname)
		end
	end

	return pl_names, rs, Ps, i_degs, pl_refnames, st_refnames
end

# ╔═╡ af31c3ab-c459-46fd-ba5c-c0c469da5091
pl_names_orb_wakeford, rs_wakeford, Ps_wakeford, i_deg_wakeford, pl_refnames_wakeford, st_refnames_wakeford = wakeford_orb(df_ps)

# ╔═╡ 7898ea90-4863-4beb-995e-7a251235aa88
# Fixed params
df_orb = DataFrame(
	pl_name = [pl_names_orb; pl_names_orb_wakeford],
	pl_ratror = [rs; rs_wakeford],
	pl_orbper = [Ps; Ps_wakeford],
	pl_orbincl = [i_degs; i_deg_wakeford],
	pl_refnames = [pl_refnames_orb; pl_refnames_wakeford],
	st_refnames = [st_refnames_orb; st_refnames_wakeford],
)

# ╔═╡ 0f939597-807b-4381-8461-09c7b9bdf3b1
occursin_df(df, s) = df[occursin.(s, df.pl_name), :]

# ╔═╡ d6598eab-33b3-4873-b6fe-b16c6d5c37d7
max_m(p, pu, pd) = p ± mean(skipmissing((pu, abs(pd))))

# ╔═╡ 4d1a7740-24c7-4cec-b788-a386bc25f836
@mdx """
We next compute the relevant quantities for estimating atmospheric observability:
"""

# ╔═╡ 542c59fd-782f-4e15-ab6c-a450bf4714ba
compute_Mₚ(K, i, P, Mₛ) = (K/sin(i)) * cbrt(P / (2.0π*G)) * cbrt(Mₛ^2)

# ╔═╡ 7a688b2e-64bb-4bc5-a799-4b267e5c30ad
compute_aRₛ(ρₛ, P) = cbrt((G * P^2 * ρₛ)/(3.0π))

# ╔═╡ 29603a24-316b-4d01-9605-6d49424fc7ff
compute_aRₛ(2.1070418u"g/cm^3", 0.813475u"d") |> NoUnits

# ╔═╡ ee3990d1-91b6-4c47-bba1-d016c75476da
compute_aRₛ(1.4945u"Msun/Rsun^3", 0.813475u"d") |> u"Msun/Msun"

# ╔═╡ 763503a1-9c3b-4353-b396-96e62c48c2be
compute_Tₚ(Tₛ, aRₛ; α=0.0) = Tₛ * (1.0 - α)^(1//4) * (0.5/aRₛ)^(1//2)

# ╔═╡ 47831596-0483-4420-a071-832183b1c3bb
@mdx """
## Split 'em up
"""

# ╔═╡ c7eabcc6-5139-448d-abdb-ec752788bd59
strip_u(u) = x -> ustrip(u, x)

# ╔═╡ e0365154-d6c8-4db2-bb85-bf2536a3aa74
function compute_Teq(T, R, a; α)
	T * (1.0 - α)^0.25 * sqrt(0.5 * R / a)
end

# ╔═╡ 05d65745-6972-41fe-8308-e5c97c85692b
get_gₚ(Mₚ, RₚRₛ, Rₛ) = G * Mₚ / (RₚRₛ^2 * Rₛ^2)

# ╔═╡ ddd8abbb-f057-4b60-bc1b-ee7f51aaa70a
df_H2OJ = CSV.read("data/pop/H2O_J_data.csv", DataFrame;
	stripwhitespace = true,
)

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

# ╔═╡ 8d519cad-8da6-409e-b486-2bc9a6008e0f
function label_text!(ax, df, s;
	al_x=:center, al_y=:bottom, offset=nothing
)
	targ = split(s, " (")[1]
	x, y = val.(Ref(df), Ref(targ), [:pl_eqt, :pl_g])
	text!(ax, s, position=(x, y), align=(al_x, al_y), offset=offset)
end

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
# Compute TSM, assuming Mp, Rp in Earth units, Teq in K, Rₛ in solar units
function compute_TSM(Rp, Teq, Mp, Rs, J; denom=1.0)
	f = compute_scale_factor(Rp)
	return f * (Rp^3 * Teq / (Mp * Rs^2)) * 10.0^(-J/5.0) / denom
end

# ╔═╡ cff3a9c4-2a9f-4769-b309-f215bfff27f7
compute_gₚ(Mₚ, RₚRₛ, Rₛ) = G * Mₚ / (RₚRₛ^2 * Rₛ^2)

# ╔═╡ 1e587a84-ed43-4fac-81cf-134a4f3d65d5
compute_H(μ, T, g) = k * T / (μ * g)

# ╔═╡ c1e63cf3-7f30-4858-bdd6-125d2a99529f
compute_ΔD(H, Rₚ, Rₛ) = 2.0 * H * Rₚ/Rₛ^2

# ╔═╡ 463341ab-318d-402f-9545-b44cf19a75ea
df_complete = let
	df = leftjoin(leftjoin(df_K, df_orb; on=:pl_name), df_exoarchive_TICv8;
		on = :pl_name
	) |> dropmissing!
	@transform! df begin
		:pl_massj = @. compute_Mₚ(
			:pl_rvamp, :pl_orbincl, :pl_orbper, :st_mass
		) |> u"Mjup"
		:pl_radj = @. (:pl_ratror * :st_rad) |> u"Rjup"
		:pl_ratdor = @. compute_aRₛ(:st_dens_solar, :pl_orbper) |> NoUnits
	end
	@transform! df begin
		:pl_g_SI = @. compute_gₚ(
			:pl_massj, :pl_ratror, :st_rad
		) |> u"m/s^2"
		:pl_eqt = @. compute_Tₚ(:st_teff, :pl_ratdor)
	end
	@transform! df begin
		:pl_H_km = @. compute_H(
			2.0*u"u", :pl_eqt, :pl_g_SI
		) |> u"km"

		# Assumes Mp, Rp in Earth units, Teq in K, Rₛ in solar units
		:TSM = @. compute_TSM(
			:pl_radj |> strip_u(u"Rearth"),
			:pl_eqt |> strip_u(u"K"),
			:pl_massj |> strip_u(u"Mearth"),
			:st_rad |> strip_u(u"Rsun"),
			:sy_jmag
		)
	end
	@chain df begin
		@aside begin
			targ_idx = _.pl_name .== "HAT-P-23 b"
			TSM_target = _.TSM[targ_idx][1]
		end
		@transform! begin
			:TSMR = :TSM / TSM_target
	
			:ΔD_ppm = @. compute_ΔD(
				:pl_H_km, :pl_radj, :st_rad
			) * 5.0 * 1e6 |> NoUnits
		end
		@aside ΔD_ppm_targ = _.ΔD_ppm[targ_idx][1]
		@transform! begin
			:ΔDR_ppm = :ΔD_ppm / ΔD_ppm_targ
		end
	end
	# @select!(df,
	# 	:TIC, :pl_name,
	# 	:Rₛ, :T_eff, :ρₛ, :Mₛ,
	# 	:K_mps,
	# 	:r, :aRₛ, :P_d, :i_deg,
	# 	:Mₚ_J, :Rₚ_J, :Tₚ, :gₚ_SI,
	# )
end

# ╔═╡ 9d7c6dff-9662-4752-8b0f-cd30357d2d85
occursin_df(df_complete, "WASP-43")[:, [:st_rad, :pl_eqt, :pl_massj]]

# ╔═╡ 373e3a8c-39f8-4656-9fcb-e0fc21cce353
df_all = @chain df_complete begin
	@select begin
		:pl_name
		:pl_eqt = @. value(:pl_eqt) |> strip_u(u"K")
		:pl_eqt_err = @. uncertainty(:pl_eqt) |> strip_u(u"K")
		:pl_g = @. value(:pl_g_SI) |> strip_u(u"m/s^2")
		:pl_g_err = @. uncertainty(:pl_g_SI) |> strip_u(u"m/s^2")
		:TSMR = value.(:TSMR)
	end
end

# ╔═╡ 998af70c-d784-4791-9261-a6dcbec8c824
df_HGHJ = @chain df_complete begin
	@rsubset (1.0 ≤ :TSMR) & (20.0u"m/s^2" ≤ :pl_g_SI) & (0.89u"Rjup" ≤ :pl_radj)
	sort(:TSMR, rev=true) # To stack smaller circles on top in Figure
	@select begin
		:pl_name
		:pl_eqt = @. value(:pl_eqt) |> strip_u(u"K")
		:pl_eqt_err = @. uncertainty(:pl_eqt) |> strip_u(u"K")
		:pl_g = @. value(:pl_g_SI) |> strip_u(u"m/s^2")
		:pl_g_err = @. uncertainty(:pl_g_SI) |> strip_u(u"m/s^2")
		:TSMR = value.(:TSMR)
		:pl_refnames
		:st_refnames
	end
end

# ╔═╡ 157a44f4-6191-4407-98a4-3c8c43817a65
df_HGHJ_no_H2OJ = filter(x -> x.pl_name ∈ ["HAT-P-23 b", "WASP-50 b"], df_HGHJ)

# ╔═╡ 94a5f868-d043-4c1f-831c-17ebabd3df6c
latextabular(df_HGHJ, latex=false, fmt="%.2f") |> PlutoUI.Text

# ╔═╡ b3ae27e9-2564-4f4c-8c51-5a40b2705ecf
df_HGHJ_paper = @chain df_complete begin
	@rsubset (1.0 ≤ :TSMR) & (20.0u"m/s^2" ≤ :pl_g_SI) & (0.89u"Rjup" ≤ :pl_radj)
	sort(:TSMR, rev=true) # To stack smaller circles on top in Figure
	@select begin
		:pl_name
		:pl_eqt
		:pl_g = @. value(:pl_g_SI) |> strip_u(u"m/s^2")
		:TSMR = value.(:TSMR)
		:pl_refnames
		:st_refnames
	end
end

# ╔═╡ d6449d05-ee95-4bda-8636-37c71e422944
df_wakeford = let
	df = leftjoin(df_H2OJ, df_complete, on=:pl_name)
	@select df begin
		:pl_name = :pl_name
		:pl_eqt = @. value(:pl_eqt) |> strip_u(u"K")
		:pl_eqt_err = @. uncertainty(:pl_eqt) |> strip_u(u"K")
		:pl_g = @. value(:pl_g_SI) |> strip_u(u"m/s^2")
		:pl_g_err = @. uncertainty(:pl_g_SI) |> strip_u(u"m/s^2")
		:H2OJ
		:H2OJ_unc
		:TSMR = value.(:TSMR)
	end
end

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
				textsize = 20,
				font = AlgebraOfGraphics.firasans("Medium"),
			),
			Lines = (linewidth=3,),
			Scatter = (linewidth=10,),
			Text = (font = AlgebraOfGraphics.firasans("Regular"), textsize=20),
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
		:pl_g => "Surface gravity (m/s²)",
	) *
	(
		data(df_all) * visual(color=(:darkgrey, 0.25))
		
		+ data(df_wakeford)
			* mapping(color=:H2OJ => "H₂O - J")
			* visual(
				marker=:rect, markersize=20, strokewidth=1, colormap=:cividis
			) #+
		
		+ data(df_HGHJ_no_H2OJ) * visual(marker='□', markersize=20)
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
	
	label_text!(ax, df_wakeford, "WASP-43 b (Weaver+ 2020)";
		al_x=:left, offset=(10, 0)
	)

	label_text!(ax, df_wakeford, "HD 189733 b (Sing+ 2016)";
		al_x=:left, offset=(0, 8)
	)
	label_text!(ax, df_HGHJ_no_H2OJ, "HAT-P-23 b (Weaver+ 2021)";
		al_x=:left, offset=(0, 8)
	)
	label_text!(ax, df_HGHJ_no_H2OJ, "WASP-50 b (this work)";
		al_x=:right, offset=(8, 8)
	)
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
		:pl_g => "Surface gravity (m/s²)",
		#color = :ΔD_ppm => "ΔD (ppm)",
		markersize = :TSMR => (x -> markersize_factor*x),
	)
	p = m * data(df_HGHJ) * visual(colormap=:viridis, marker='○')

	fg = draw(p;
		axis = (; limits=((0, 3_400), (-1, 55)), yticks=0:10:50),
		figure = (; resolution=FIG_LARGE),
	)
	ax = fg.grid[1].axis

	# HGHJ g boundary
	hl = hlines!(ax, 20.0, color=:darkgrey, linestyle=:dash)
	translate!(hl, 0, 0, -1) # Place behind plot markers

	# Annotate HGHJs with tspec observations
	HP23x, HP23y = val.(Ref(df_HGHJ), Ref("HAT-P-23 b"), [:pl_eqt, :pl_g])
	annotate_text!(
		ax,
		"HAT-P-23b",
		(HP23x[1], HP23y[1]) .- (-100, -5),
		(HP23x[1], HP23y[1]),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	W43x, W43y = val.(Ref(df_HGHJ), Ref("WASP-43 b"), [:pl_eqt, :pl_g])
	annotate_text!(
		ax,
		"WASP-43b",
		(W43x[1], W43y[1]) .- (300, -2),
		(W43x[1], W43y[1]),
		0.0,
		0.0;
		align = (:center, :top),
	)
	W50x, W50y = val.(Ref(df_HGHJ), Ref("WASP-50 b"), [:pl_eqt, :pl_g])
	annotate_text!(
		ax,
		"WASP-50b",
		(W50x[1], W50y[1]) .- (200, -2.5),
		(W50x[1], W50y[1]),
		0.5,
		0.1;
		align = (:center, :baseline),
	)
	hdx, hdy = val.(Ref(df_HGHJ), Ref("HD 189733 b"), [:pl_eqt, :pl_g])
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
# ╠═fccd3615-7726-45ac-8d5d-c4f1fc4d5425
# ╠═0d629db3-7370-406f-989b-7a2caca020dc
# ╠═380d05a4-35e9-4db4-b35a-b03de9e695ee
# ╠═f8b4def8-a46f-4cbc-83c0-ff44a39c1571
# ╠═86a99042-bb9b-43e6-87ae-d76f88b10533
# ╠═92fbb7d7-9782-44d4-b1b7-6db89d78a032
# ╟─97c9f1ae-21da-4f99-94ff-a8adaabf30bb
# ╠═2584860a-8e24-49f7-a7d5-4c99c8deda8e
# ╠═759b0ca7-ade4-4929-afa5-51e0ab133a5b
# ╠═2ae250e2-cc4d-4824-9795-a8bc0a4b469b
# ╟─3ed05d01-b489-46d6-bcd4-9028d756ab35
# ╠═0aa8aaf2-5343-4b6e-a47b-cf0fc8d27643
# ╠═0c793036-d7c6-4a56-9ef6-f58b02e6530c
# ╠═2d63caf3-dd64-483d-8c85-5085d7aad2ac
# ╠═f306bc5a-5597-4a52-b6e3-364a4230367d
# ╠═af31c3ab-c459-46fd-ba5c-c0c469da5091
# ╠═7898ea90-4863-4beb-995e-7a251235aa88
# ╠═463341ab-318d-402f-9545-b44cf19a75ea
# ╠═9d7c6dff-9662-4752-8b0f-cd30357d2d85
# ╠═0f939597-807b-4381-8461-09c7b9bdf3b1
# ╠═d6598eab-33b3-4873-b6fe-b16c6d5c37d7
# ╟─4d1a7740-24c7-4cec-b788-a386bc25f836
# ╠═542c59fd-782f-4e15-ab6c-a450bf4714ba
# ╠═7a688b2e-64bb-4bc5-a799-4b267e5c30ad
# ╠═29603a24-316b-4d01-9605-6d49424fc7ff
# ╠═ee3990d1-91b6-4c47-bba1-d016c75476da
# ╠═763503a1-9c3b-4353-b396-96e62c48c2be
# ╠═47831596-0483-4420-a071-832183b1c3bb
# ╠═373e3a8c-39f8-4656-9fcb-e0fc21cce353
# ╠═998af70c-d784-4791-9261-a6dcbec8c824
# ╠═b3ae27e9-2564-4f4c-8c51-5a40b2705ecf
# ╠═4e6c0113-255f-46fb-84ad-f8f30b79a219
# ╠═d6449d05-ee95-4bda-8636-37c71e422944
# ╠═c7eabcc6-5139-448d-abdb-ec752788bd59
# ╠═e0365154-d6c8-4db2-bb85-bf2536a3aa74
# ╠═05d65745-6972-41fe-8308-e5c97c85692b
# ╠═ddd8abbb-f057-4b60-bc1b-ee7f51aaa70a
# ╠═157a44f4-6191-4407-98a4-3c8c43817a65
# ╠═c0f576a7-908d-4f10-86e7-cadbb7c77c09
# ╠═8d519cad-8da6-409e-b486-2bc9a6008e0f
# ╠═c1cd9292-28b9-4206-b128-608aaf30ff9c
# ╟─0f9262ef-b774-45bc-bdab-46860779683d
# ╠═94a5f868-d043-4c1f-831c-17ebabd3df6c
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
