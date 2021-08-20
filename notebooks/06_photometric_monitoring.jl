### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 9e2ce576-c9bd-11eb-0699-47af13e79589
begin
	import Pkg
	Pkg.activate(Base.current_project())
	
	using AlgebraOfGraphics
	using CSV
	using CairoMakie
	using CCDReduction: fitscollection
	using Colors
	using DataFrames
	using DataFramesMeta
	using Dates
	using DelimitedFiles
	using Glob
	using ImageFiltering
	using KernelDensity
	using Latexify
	using LaTeXStrings
	using LombScargle
	using Measurements
	using Measurements: value, uncertainty
	using NaturalSort
	using OrderedCollections
	using Printf
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
	using Unitful
	
	# Python setup
	#ENV["PYTHON"] = "~/miniconda3/envs/WASP-50b/bin/python"
	Pkg.build("PyCall")
	using PyCall
	
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	
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

# ╔═╡ 670b88e4-1e96-494d-bfcc-235092bb6e96
md"""
# Photometric Monitoring

In this notebook we gather and analyze the available photometric data for this target.

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 0cbe4263-799f-4ee3-9a94-3ba879528b01
md"""
## ASAS-SN 🌍

We start by downloading [all photometric monitoring data](https://asas-sn.osu.edu/photometry/7df0eb29-0b68-57ef-8ce2-83dc7b5684da) gathered by the [ASAS-SN](https://asas-sn.osu.edu/) survey, which we include below.
"""

# ╔═╡ b00c28a2-26b1-442e-a347-39fb66b825a0
md"""
### Data inspection

Let's start by checking the data for any corrupted or missing data points:
"""

# ╔═╡ fa233a2c-7e89-4e71-85b8-824c5c341650
df_ASASSN = CSV.File(
	"data/photometric/ASAS-SN/AP37847073.csv",
	normalizenames = true,
) |> DataFrame

# ╔═╡ 3345ff07-b04c-498e-8161-f87082448137
nrow(df_ASASSN)

# ╔═╡ 9094c4a4-3b75-4e21-97a7-600de734867b
describe(df_ASASSN, :all)

# ╔═╡ 2b2dd83b-ce99-4551-b8aa-79ca6db5dd06
function name(dt_julian)
	dt = julian2datetime(dt_julian)
	return "$(year(dt)) $(monthname(dt)) $(day(dt))"
end

# ╔═╡ 6eaf882c-0cb5-415f-b8fe-c071ee25a895
md"""
Looks good. According to the table above, the data spans from
**$(join(name.(extrema(df_ASASSN.hjd)), " - "))** from two cameras (**bd**, commissioned 2013 December at Haleakala; **bh**, commissioned 2015 July at CTIO), both in the V band:
"""

# ╔═╡ 98704543-8cb7-4fca-b601-2a8d2dfa4833
begin
	# Define transit epochs
	utc_transit_dates = ["2013-12-19", "2015-09-27", "2016-12-11"]
	transit_dates = DateTime.(utc_transit_dates)
	julian_transit_dates = datetime2julian.(transit_dates) |> collect
end;

# ╔═╡ 8e6008ca-a762-4450-a09e-bd0c2bbac4f2
let
	phot_mon = data(df_ASASSN) * mapping(
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
			position = Point2f0(jd, 11.83),
			textsize = 14,
			align = (:left, :center),
			offset = Point2f0(10, 0),
			color = :grey,
		)
	end
		
	#save("phot_mon.png", fig, px_per_unit=3)
	
	fig #|> as_svg
end

# ╔═╡ 7e5ea835-8eb2-44d5-905d-433767b6b91a
md"""
Since these observations share the same filter, we will treat them as if they came from a single instrument. Next, we will convert the recorded HJD times to BTJD (BJD - 2457000), to place them on the same time axis as the TESS data.
"""

# ╔═╡ 5bc41820-3782-4f82-80b6-6da62805ca8f
md"""
### Time conversion

We will make use of astropy's coordinates and time modules to first make the conversion from HJD to BJD:
"""

# ╔═╡ 5bf1136b-2d13-4463-8d74-9ade1e2cee96
# begin
# 	py"""
# 	from astropy.coordinates import SkyCoord, EarthLocation
# 	from astropy import units as u
# 	from astropy.time import Time
	
	
# 	# https://gist.github.com/StuartLittlefair/4ab7bb8cf21862e250be8cb25f72bb7a
# 	def helio_to_bary(coords, hjd, obs_name):
# 	    helio = Time(hjd, scale="utc", format="jd")
# 	    obs = EarthLocation.of_site(obs_name)
# 	    star = SkyCoord(coords, unit=(u.hour, u.deg)) 
# 	    ltt = helio.light_travel_time(star, "heliocentric", location=obs)
# 	    guess = helio - ltt
	    
# 	    # If we assume guess is correct - how far is heliocentric time away
# 	    # from true value?
# 	    delta = (
# 	    guess + guess.light_travel_time(star, "heliocentric", obs)
# 	    ).jd - helio.jd
	    
# 	    # Apply this correction
# 	    guess -= delta * u.d
# 	    ltt = guess.light_travel_time(star, 'barycentric', obs)
	    
# 	    return guess.tdb + ltt
# 	"""
# 	helio_to_bary(coords, hjd, obs_name) = py"helio_to_bary"(
# 		coords, hjd, obs_name
# 	).value
# end

# ╔═╡ c86f5adc-6e17-44c4-b754-1b5c42557809
const coords_W50 = "1:12:43.2 +31:12:43" # RA, Dec of WASP-50

# ╔═╡ 5f1ca89d-61e8-4ade-93e1-13716dc5fd46
loc_to_BJD(t, coords, camera) = camera == "Haleakala" ?
	helio_to_bary(coords, t, "Haleakala") :
	helio_to_bary(coords, t, "CTIO")

# ╔═╡ 42b72605-f89e-42c4-a159-d43e4620140f
# df_ASASSN[!, :bjd] = loc_to_BJD.(df_ASASSN.hjd, coords_W50, df_ASASSN.camera);

# ╔═╡ 79d08932-66f3-4ed9-bc13-f1ac3229e95d
df_ASASSN

# ╔═╡ 92548af3-9a26-4202-88f2-ba3a31181686
# begin
# 	df_sorted = sort(df_ASASSN, :bjd)
# 	t_ASASSN, f_ASASSN, f_err_ASASSN = eachcol(
# 		df_sorted[!, [:bjd, :flux_mJy_, :flux_err]]
# 	)
# 	#t_ASASSN .-= 2.457e6
# end;

# ╔═╡ e36d1322-5aa4-4513-bb89-410a4bb6b750
# t_ASASSN_binned, f_ASASSN_binned, f_ASASSN_binned_err = bin_lc(
# 	t_ASASSN, f_ASASSN, f_err_ASASSN, t_window
# );

# ╔═╡ 7f864c2d-e6a2-4774-b56e-6131041c3a00
# df_sorted |> describe

# ╔═╡ dbe317fe-540d-44e8-b8e7-6c465c79559f
md"""
``t_\text{window}`` = $(@bind t_window Slider(1:30, default=7, show_value=true)) days
"""

# ╔═╡ 9a195365-e68f-43e2-8870-c09153e2ce91
md"""
### Plot

With the BJD times computed, we can now plot the ASAS-SN photometry, binned to **$(t_window) days**:
"""

# ╔═╡ f3425d9c-861e-4b26-b352-bd0669c7f1f9
# let
# 	fig = Figure(resolution=(800, 800))
	
# 	### Photometry plot ####
# 	ax = Axis(fig[1, 1], xlabel="Time (BTJD)", ylabel="Relative flux (ppm)")
	
# 	# Mark transit epochs
# 	Δjulian_transit_dates = julian_transit_dates .- 2.457e6
# 	vlines!(ax, Δjulian_transit_dates;
# 		linestyle = :dash,
# 		color = :darkgrey,
# 	)
	
# 	# Label transit epochs
# 	for (i, (utc, jd)) in enumerate(zip(utc_transit_dates, Δjulian_transit_dates))
# 		text!(ax, "Transit $i\n$utc";
# 			position = Point2f0(jd, 5.0e4),
# 			textsize = 14,
# 			align = (:left, :center),
# 			offset = Point2f0(10, 0),
# 			color = :grey,
# 		)
# 	end
	
# 	ax_phot, t_binned, f_binned, f_err_binned, lc_binned, f_med = plot_phot!(
# 		ax, t_ASASSN, f_ASASSN, f_err_ASASSN;
# 		t_offset=2.457e6, relative_flux=true, binsize=t_window
# 	)
	
# 	#### Periodogram #######
# 	ax_pg = Axis(fig[2, 1], xlabel="Periods (days)", ylabel="log10 Power")
# 	Ps, powers, faps = compute_pg(lc_binned, 5, 30)
# 	lines!(ax_pg, Ps, log10.(powers*f_med/1e6))
# 	hlines!(ax_pg, log10.(faps), color=:darkgrey, linestyle=:dash, label="\n\n[1, 5, 10]% FAPs")
# 	axislegend()
	
# 	CSV.write(
# 		"/home/mango/Desktop/WASP50LC_ASASSN_binned.csv",
# 		DataFrame(:t=>t_binned, :f=>f_binned, :f_err=>f_err_binned),
# 	)
			
# 	fig
# end

# ╔═╡ 682499cb-af79-48f7-9e74-0185894f65fe
#= md"""
For comparison, the cadence for this data has an average of $(t_ASASSN |> diff |> mean) days.
""" =#

# ╔═╡ 78d85c7c-da06-41ab-915b-48d93a010967
md"""
## TESS 🌌
We next turn to the TESS photometry.
"""

# ╔═╡ 97e7feee-11b2-4a35-9327-b5c0d05b2a23
md"""
### Light curve collection

First we use [`lightkurve`](https://docs.lightkurve.org/whats-new-v2.html) to download the available data products from TESS:
"""

# ╔═╡ 0c790d2f-64d4-4e13-9629-a9725cd7086d
lk = pyimport("lightkurve")

# ╔═╡ 708d54a5-95fd-4f15-9681-f6d8e7b9b05c
all_srs = lk.search_lightcurve("WASP-50")

# ╔═╡ dff46359-7aec-4fa1-bc7a-89785dfca0e8
srs = lk.search_lightcurve("WASP-50", author=["SPOC"], exptime=120)

# ╔═╡ 34fcd73d-a49c-4597-8e63-cfe2495eee48
md"""
From the $(length(all_srs)) data products found, we see that $(length(srs)) are available from the [Science Processing Operations Center (SPOC)](https://heasarc.gsfc.nasa.gov/docs/tess/pipeline.html):
"""

# ╔═╡ 6e62b3cc-96ce-43fd-811b-4b2b102cfd61
lcs = srs.download_all(flux_column="pdcsap_flux")

# ╔═╡ 241c462c-3cd9-402d-b948-b9b1f608b727
md"""
We show the normalized PDCSAP flux below for each sector: 
"""

# ╔═╡ ec12acb8-9124-4cc0-8c9f-6525c1565dfd
begin
	py"""
	def oot_flux(lc, P, t_0, dur):
		in_transit = lc.create_transit_mask(P, t_0, dur)
		lc_oot = lc[~in_transit]
		return lc_oot
	"""
	oot_flux(lc, P, t_0, dur) = py"oot_flux"(lc, P, t_0, dur)
end

# ╔═╡ 2952e971-bce5-4a1e-98eb-cb2d45c8c5a8
Time = pyimport("astropy.time").Time

# ╔═╡ 31d5bc92-a1f2-4c82-82f2-67755f9aa235
begin
	P = 1.9550931258
	t_0 = Time(2455558.61237, format="jd")
	dur = 1.83 * (1.0 / 24.0)
	
	lcs_cleaned = []
	lcs_oot = []
	for lc in lcs
		lc_cleaned = lc.remove_nans().normalize()
		lc_oot = oot_flux(lc_cleaned, P, t_0, dur).remove_outliers(sigma=3.0)
		push!(lcs_cleaned, lc_cleaned)
		push!(lcs_oot, lc_oot)
	end
	push!(lcs_cleaned, lk.LightCurveCollection([lcs_cleaned...]).stitch())
	push!(lcs_oot, lk.LightCurveCollection([lcs_oot...]).stitch())
end;

# ╔═╡ 82222ee8-f759-499d-a072-c219cc33ccad
let
	fig = Figure(resolution=FIG_TALL)
	
	for (i, (lc, lc_oot)) ∈ enumerate(zip(lcs_cleaned, lcs_oot))
		ax = Axis(fig[i, 1])
		errorbars!(
			ax,
			lc.time.value,
			lc.flux,
			lc.flux_err,
			label = """
			Sector $(lc.meta["SECTOR"]), $(lc.meta["AUTHOR"])
			"""
		)
		#scatter!(fig[i, 1], lc.time.value, lc.flux)
		
		errorbars!(ax, lc_oot.time.value, lc_oot.flux, lc_oot.flux_err;
			color = COLORS[2], label="oot",
			markersize = 5,
		)

		axislegend()
	end

	#linkyaxes!(filter(x -> x isa Axis, fig.content)...)

	Label(fig[end+1, 1], "Time (BTJD days)", tellwidth=false)
	Label(fig[1:end-1, 0], "Relative flux", rotation=π/2)
	
	fig
end

# ╔═╡ 897da774-0e18-43f8-8337-4048f12d8e45
# oot_flux(lcs[2], P, dur, t_0).to_csv("/home/mango/Desktop/WASP50_S31_oot.csv")

# ╔═╡ e343fcb5-f633-46eb-a534-7772e11c3f60
# oot_flux(lcs.stitch(), P, dur, t_0).to_csv("/home/mango/Desktop/WASP50_combined_oot.csv")

# ╔═╡ bf01b1c2-e4c6-4bed-8e79-a20f273cc387
function plot_photometry(lcs, yfield)
	fig = Figure(resolution=FIG_TALL)
	
	for (i, lc) in enumerate(lcs)
		lines(
			fig[i, 1],
			lc.time.value,
			getproperty(lc, yfield).value,
			label = """
			Sector $(lc.meta["SECTOR"]), $(lc.meta["AUTHOR"])
			$(lc.meta["EXPOSURE"]) s
			$(srs[i].exptime)
			"""
		)
		
		axislegend()
	end
	
	linkyaxes!(filter(x -> x isa Axis, fig.content)...)
			
	if yfield == :sap_flux
		ylabel = "SAP Flux (e⁻/s)"
	elseif yfield == :pdcsap_flux
		ylabel = "PDCSAP Flux (e⁻/s)"
	else
		ylabel = "Flux"
	end
	
	Label(fig[end+1, 1], "Time (BTJD days)", tellwidth=false)
	Label(fig[1:end-1, 0], ylabel, rotation=π/2)
	
	fig
end

# ╔═╡ 99fcf959-665b-44cf-9b5f-fd68a919f846
md"""
### Periodogram
"""

# ╔═╡ a203f10c-7b7f-4b2f-b020-4154138ce5e5
md"""
#### `tess_sip`
"""

# ╔═╡ d1f7ed4b-4599-48bd-aac5-93920dae9151
# SIP = pyimport("tess_sip").SIP
# binsize = 1.2 #12 * 3600 * 1 / 86_400
# lcs_combined = lcs #lk.LightCurveCollection([lc.bin(binsize) for lc ∈ lcs])
# lcs_S04 = lk.LightCurveCollection([lcs_combined[1]])
# lcs_S31 = lk.LightCurveCollection([lcs_combined[2]])
# SIP_kwargs = (min_period=1.0, max_period=35.0, nperiods=100)
# r_S04, r_S31, r_combined = SIP.((lcs_S04, lcs_S31, lcs_combined); SIP_kwargs...)

# function plot_SIP_flux(r)
# 	lc_raw = r["raw_lc"]
# 	lc_corr = r["corr_lc"]
	
# 	fig = Figure()
# 	ax = Axis(fig[1, 1])

# 	lines!(ax, lc_raw.time.value, lc_raw.flux, color=:darkgrey)
# 	scatter!(ax, lc_corr.time.value, lc_corr.flux)

# 	#xlims!(ax, 1400, 1450)
# 	ylims!(ax, 0.95, 1.05)
	
# 	fig
# end

# pgram_S04, plan_S04 = compute_pgram(r_S04["corr_lc"])
# pgram_S31, plan_S31 = compute_pgram(r_S31["corr_lc"])
# pgram_combined, plan_combined = compute_pgram(r_combined["corr_lc"])

# ╔═╡ 5a275d11-f503-4623-b200-097abdfb2a76
lcs_oot[1].bin(0.5)

# ╔═╡ de104bdf-e95a-4a6b-9178-c6b2a2f2f5ea
function compute_pgram(lc; min_period=0.5, max_period=30)
	plan = LombScargle.plan(
		lc.time.value, lc.flux .± lc.flux_err,
		minimum_frequency = 1.0 / max_period,
		maximum_frequency = 1.0 / min_period,
	)
	return lombscargle(plan), plan
end

# ╔═╡ d7f034c5-5925-4b91-9bea-1068a7ce9252
begin
	pgrams, plans, P_maxs = [], [], []
	for lc in lcs_oot
		pgram, plan = compute_pgram(lc)
		P_max = findmaxperiod(pgram)[1]
		
		push!(pgrams, pgram)
		push!(plans, plan)
		push!(P_maxs, P_max)
	end
end

# ╔═╡ bdee2a77-5554-412e-b461-1f548e9d880b
compute_pgram(lcs_oot[1].bin(0.005))

# ╔═╡ 94d05a5b-b05e-4407-bcd3-7d625680a262
let
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="Period (days)", ylabel="Normalized power")
	
	for (pgram, plan, P_max) in zip(pgrams, plans, P_maxs)
		# Compute FAPs
		b = LombScargle.bootstrap(100, plan)
		
		# Plot
		lines!(ax, periodpower(pgram)..., label="$(round(P_max, digits=2))")
		hlines!(ax, collect(fapinv.(Ref(b), (0.01, 0.05, 0.1))))
	end
	
	axislegend()
	
	fig
end

# ╔═╡ a50ef756-ade6-48a3-8d3a-17b56ce03c26
md"""
### Folded lightcurve
"""

# ╔═╡ 3128e57f-df4f-4811-b867-8a293d7d536d
function compute_pgram_model(lc, P)
	t_fit = lc.time.value
	s_fit = LombScargle.model(
		lc.time.value,
		lc.flux,
		lc.flux_err,
		1 / P,
	)
	lc_fit = lk.LightCurve(t_fit, s_fit)
	lc_fit_folded = lc_fit.fold(P)
	return lc_fit_folded
end

# ╔═╡ 97ced6ba-ff74-46b4-90d5-18e7b2f1b903
begin
	lcs_folded = []
	lcs_folded_binned = []
	lcs_fit_folded = []
	for (lc, P) in zip(lcs_oot, P_maxs)
		# Data
		lc_folded = lc.fold(P)
		push!(lcs_folded, lc_folded)
		push!(lcs_folded_binned, lc_folded.bin(P))
		
		# Model
		lc_fit_folded = compute_pgram_model(lc, P)	
		push!(lcs_fit_folded, lc_fit_folded)
	end
end

# ╔═╡ 49bcddbe-d413-48ae-91d8-92bcebf40518
let
	fig = Figure()
	
	axs = []
	for (i, (lc_folded, lc_folded_binned, lc_fit_folded)) in enumerate(zip(
				lcs_folded, lcs_folded_binned, lcs_fit_folded
		))
		ax = Axis(fig[i, 1])
		push!(axs, ax)
		scatter!(ax, lc_folded.time.value, lc_folded.flux)
		scatter!(ax, lc_folded_binned.time.value, lc_folded_binned.flux)
		lines!(ax, lc_fit_folded.time.value, lc_fit_folded.flux, color=:lightgreen)
	end
	
	linkaxes!(axs...)
	
	fig
end

# ╔═╡ 87e7f166-f022-4b53-a3a3-c96b9d377781
# let
# 	tmin, tmax = extrema(lc_corr.time.value)
# 	Δt = median(diff(lc_corr.time.value))
# 	t = tmin:0.001:tmax
# 	s = oneunit.(t)
# 	s_err = median(lc_corr.flux_err) .* s
	
# 	pgram = lombscargle(t, s .± s_err)
	
# 	lines(periodpower(pgram)...)
# end

# ╔═╡ 056281a2-4786-45eb-a9fa-57515153f66c
md"""
#### Spot parameter estimation
"""

# ╔═╡ 3a612743-7071-4d85-a48d-0a4b12facffc
ΔL = minimum(lc_folded.flux) / median(lc_folded.flux)

# ╔═╡ 278ea804-e2dd-4ca8-9a20-0e9a25746e02
T₀ = 5_520

# ╔═╡ d26122f1-1602-440c-8ba9-72469a782104
f_sp(T_sp, T₀, ΔL) = (T₀^4 /  (T₀^4 - T_sp^4)) * (1 - ΔL)

# ╔═╡ c0d53fc0-467b-4e49-a829-f149e14e3d08
f_sp.((2_200, 2_800), T₀, ΔL) .* 100

# ╔═╡ a62cae71-f73f-49bc-992c-ba7dbf4792d9
md"""
With the TESS data in hand, we now run it through [TESS_SIP](https://github.com/christinahedges/TESS-SIP) to analyze the final periodogram:
"""

# ╔═╡ 1a86efc8-0bb6-44e4-8568-82fdb3409c25
md"""
### Window function
"""

# ╔═╡ 2215ed86-fa78-4811-88ab-e3521e4a1dea
function compute_window_func(lc; P_min=5, P_max=30)
	t = lc.time.value
	t_start, t_end = t |> extrema
	Δt = median(diff(t))
	t_uniform = t #t_start:Δt:t_end
	
	f = oneunit.(t_uniform)
	f_err = median(lc.flux_err) .* f
	lc_window_func = lk.LightCurve(time=t_uniform, flux=f, flux_err=f_err)
	
	return compute_pgram(lc_window_func) # period, power, faps
end

# ╔═╡ a675c241-dba2-4764-8e37-ca24ab98f8cd
ppgram, pplan = compute_window_func(lc_corr)

# ╔═╡ 4377b1c0-bd0e-4f48-ad34-4aa4ae297aab
let
	fig = Figure(); ax = Axis(fig[1, 1])
	
	lines!(ax, periodpower(ppgram)...)
	
	#vlines!(ax, [P_max / 3.0, P_max / 2.0, P_max], color=:red)
	
	#b = LombScargle.bootstrap(100, pplan)
	#hlines!(ax, collect(fapinv.(Ref(b), (0.01, 0.05, 0.1))))
	
	fig
end

# ╔═╡ dd2d460b-1169-4760-8e74-16ea640d4448
# using Images: findlocalmaxima

# ╔═╡ ed93fd1f-79bd-4c5a-b822-680c6f6591d6
# let
# 	fig = Figure(resolution=FIG_TALL)
	
# 	P = 16.3 # Gillon+ 2011
	
# 	axs = []
# 	lcss = [lcs_S04, lcs_S31, lcs_combined]
# 	for (i, (sector, r)) in enumerate(data_dict)
# 		ax_window_func = Axis(fig[2*i-1, 1])
# 		lc = lcss[i].stitch(x -> x).remove_nans().normalize()
# 		a, b, c = compute_window_func(lc) # period, power, faps
# 		peak_idxs = findlocalmaxima(b)
# 		lines!(ax_window_func, a, b, color=:darkgrey, label="window function")
# 		scatter!(ax_window_func, a[peak_idxs], b[peak_idxs], color=COLORS[6])
		
# 		ax = Axis(fig[2*i, 1], yscale=log10)
# 		plot_periodogram!(ax, r, P)
# 		text!(ax, sector;
# 			position = (30, 4.5),
# 			textsize = 16,
# 			align = (:right, :center),
# 			color = :darkgrey,
# 		)
# 		push!(axs, ax)
# 	end
	
# 	linkaxes!(axs...)
# 	hidexdecorations!.(axs[1:2])
	
# 	Legend(fig[1:end, 2], axs[1], margin=(10, 0, 0, 0))
# 	Label(fig[end+1, 1], "Period (days)", tellwidth=false)
# 	Label(fig[1:end, 0], "log10 Power", rotation=π/2)
	
# 	fig
# end

# ╔═╡ 4fcb3ea1-bec6-4c30-82a6-4689077cdc14
# function plot_periodogram!(ax, r, P)
# 	periods = r["periods"]
#     power_source = r["power"]
#     power_background = r["power_bkg"]
#     power_relative = power_source ./ power_background

# 	lines!(ax, periods, power_relative; label="divided")
# 	lines!(ax, periods, power_source; label="source")
# 	lines!(ax, periods, power_background; label="background")
	
# 	P_max = periods[argmax(power_relative)]
# 	vlines!(ax, P, linestyle=:dash, color=:darkgrey)
# 	vlines!(ax, P/3.0, linestyle=:dash, color=:darkgrey)
# 	vlines!(ax, P_max, color=:darkgrey)
	
# 	text!(ax, "P_max = $(round(P_max, digits=2)) days";
# 		position = (30, 2),
# 		textsize = 16,
# 		align = (:right, :center),
# 		color = :darkgrey,
# 	)
# end

# ╔═╡ 18223d42-66d8-40d1-9d89-be8af46853e2
md"""
## Helper Functions
"""

# ╔═╡ 682c3732-e68f-4fdb-bd63-553223308364
begin
	Makie.convert_arguments(
		P::PointBased, v::Vector, m::AbstractVector{<:Measurement}
	) = convert_arguments(P, v, value.(m))
	
	Makie.convert_arguments(
		P::Type{<:Errorbars}, v::Vector, m::AbstractVector{<:Measurement}
	) = convert_arguments(P, v, value.(m), uncertainty.(m))	
end

# ╔═╡ 7370a1d9-4f8e-4788-adac-b8be2bcc9643
function plot_phot!(ax, t, f, f_err; t_offset=0.0, relative_flux=false, binsize=1.0)
	t_rel = t .- t_offset
	if relative_flux
		f_med = median(f)
		Δf, Δf_err = 1e6*(f .- f_med) / f_med , 1e6*f_err / f_med
		
	else
		f_med = 1.0
		Δf, Δf_err = f, f_err
	end
	
	# Original data
	errorbars!(ax, t_rel, Δf, Δf_err, color=(:darkgrey, 0.25))
	scatter!(ax, t_rel, Δf, color=(:darkgrey, 0.25))
	
	# Binned data
	t_binned, f_binned, f_err_binned, lc_binned = bin_lc(t_rel, Δf, Δf_err, binsize)
	f_binned_err_med = median(filter(!isnan, f_err_binned))
	
	errorbars!(ax, t_binned, f_binned, f_err_binned;
		color=:grey
	)
	scatter!(ax, t_binned, f_binned;
		color=:grey, label="avg err: $(f_binned_err_med)",
	)
	
	axislegend(ax, position=:rb)
	
	return ax, t_binned, f_binned, f_err_binned, lc_binned, f_med
end

# ╔═╡ ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
md"""
## Notebook setup
"""

# ╔═╡ 4475d6c8-c6cc-4a5a-903e-96aa0c457864
let
	t = range(0.01, stop = 10pi, length = 1000) # Observation times
	s = sinpi.(t) .+ 1.2cospi.(t) .+ 0.3rand(length(t)) # The noisy signal
	
	# Pick-up the best frequency
	f = findmaxfreq(lombscargle(t, s, maximum_frequency=10, samples_per_peak=20))[1]
	t_fit = range(0, stop = 1, length = 50)
	s_fit = LombScargle.model(t, s, f, t_fit/f) # Determine the model
	
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	scatter!(ax, mod.(t.*f, 1), s, label="Phased data", title="Best Lomb-Scargle frequency: $f")
	lines!(t_fit, s_fit, label="Best-fitting model", linewidth=4, color=:orange)
	
	axislegend()
	
	fig
end

# ╔═╡ Cell order:
# ╟─670b88e4-1e96-494d-bfcc-235092bb6e96
# ╟─0cbe4263-799f-4ee3-9a94-3ba879528b01
# ╟─b00c28a2-26b1-442e-a347-39fb66b825a0
# ╠═fa233a2c-7e89-4e71-85b8-824c5c341650
# ╠═3345ff07-b04c-498e-8161-f87082448137
# ╠═9094c4a4-3b75-4e21-97a7-600de734867b
# ╟─6eaf882c-0cb5-415f-b8fe-c071ee25a895
# ╟─2b2dd83b-ce99-4551-b8aa-79ca6db5dd06
# ╠═98704543-8cb7-4fca-b601-2a8d2dfa4833
# ╠═8e6008ca-a762-4450-a09e-bd0c2bbac4f2
# ╟─7e5ea835-8eb2-44d5-905d-433767b6b91a
# ╟─5bc41820-3782-4f82-80b6-6da62805ca8f
# ╠═5bf1136b-2d13-4463-8d74-9ade1e2cee96
# ╠═c86f5adc-6e17-44c4-b754-1b5c42557809
# ╠═5f1ca89d-61e8-4ade-93e1-13716dc5fd46
# ╠═42b72605-f89e-42c4-a159-d43e4620140f
# ╠═79d08932-66f3-4ed9-bc13-f1ac3229e95d
# ╟─9a195365-e68f-43e2-8870-c09153e2ce91
# ╠═92548af3-9a26-4202-88f2-ba3a31181686
# ╠═e36d1322-5aa4-4513-bb89-410a4bb6b750
# ╠═7f864c2d-e6a2-4774-b56e-6131041c3a00
# ╟─dbe317fe-540d-44e8-b8e7-6c465c79559f
# ╠═f3425d9c-861e-4b26-b352-bd0669c7f1f9
# ╟─682499cb-af79-48f7-9e74-0185894f65fe
# ╟─78d85c7c-da06-41ab-915b-48d93a010967
# ╟─97e7feee-11b2-4a35-9327-b5c0d05b2a23
# ╠═0c790d2f-64d4-4e13-9629-a9725cd7086d
# ╠═708d54a5-95fd-4f15-9681-f6d8e7b9b05c
# ╟─34fcd73d-a49c-4597-8e63-cfe2495eee48
# ╠═dff46359-7aec-4fa1-bc7a-89785dfca0e8
# ╠═6e62b3cc-96ce-43fd-811b-4b2b102cfd61
# ╟─241c462c-3cd9-402d-b948-b9b1f608b727
# ╠═ec12acb8-9124-4cc0-8c9f-6525c1565dfd
# ╠═2952e971-bce5-4a1e-98eb-cb2d45c8c5a8
# ╠═31d5bc92-a1f2-4c82-82f2-67755f9aa235
# ╠═82222ee8-f759-499d-a072-c219cc33ccad
# ╠═897da774-0e18-43f8-8337-4048f12d8e45
# ╠═e343fcb5-f633-46eb-a534-7772e11c3f60
# ╠═bf01b1c2-e4c6-4bed-8e79-a20f273cc387
# ╟─99fcf959-665b-44cf-9b5f-fd68a919f846
# ╟─a203f10c-7b7f-4b2f-b020-4154138ce5e5
# ╠═d1f7ed4b-4599-48bd-aac5-93920dae9151
# ╠═d7f034c5-5925-4b91-9bea-1068a7ce9252
# ╠═bdee2a77-5554-412e-b461-1f548e9d880b
# ╠═5a275d11-f503-4623-b200-097abdfb2a76
# ╠═de104bdf-e95a-4a6b-9178-c6b2a2f2f5ea
# ╠═94d05a5b-b05e-4407-bcd3-7d625680a262
# ╟─a50ef756-ade6-48a3-8d3a-17b56ce03c26
# ╠═49bcddbe-d413-48ae-91d8-92bcebf40518
# ╠═97ced6ba-ff74-46b4-90d5-18e7b2f1b903
# ╠═3128e57f-df4f-4811-b867-8a293d7d536d
# ╠═87e7f166-f022-4b53-a3a3-c96b9d377781
# ╟─056281a2-4786-45eb-a9fa-57515153f66c
# ╠═3a612743-7071-4d85-a48d-0a4b12facffc
# ╠═278ea804-e2dd-4ca8-9a20-0e9a25746e02
# ╠═d26122f1-1602-440c-8ba9-72469a782104
# ╠═c0d53fc0-467b-4e49-a829-f149e14e3d08
# ╟─a62cae71-f73f-49bc-992c-ba7dbf4792d9
# ╟─1a86efc8-0bb6-44e4-8568-82fdb3409c25
# ╠═2215ed86-fa78-4811-88ab-e3521e4a1dea
# ╠═a675c241-dba2-4764-8e37-ca24ab98f8cd
# ╠═4377b1c0-bd0e-4f48-ad34-4aa4ae297aab
# ╠═dd2d460b-1169-4760-8e74-16ea640d4448
# ╠═ed93fd1f-79bd-4c5a-b822-680c6f6591d6
# ╠═4fcb3ea1-bec6-4c30-82a6-4689077cdc14
# ╟─18223d42-66d8-40d1-9d89-be8af46853e2
# ╠═682c3732-e68f-4fdb-bd63-553223308364
# ╠═7370a1d9-4f8e-4788-adac-b8be2bcc9643
# ╟─ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
# ╠═9e2ce576-c9bd-11eb-0699-47af13e79589
# ╠═4475d6c8-c6cc-4a5a-903e-96aa0c457864
