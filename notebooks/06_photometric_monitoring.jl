### A Pluto.jl notebook ###
# v0.17.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 55beac98-0929-4a55-91f7-cee7c781498c
begin
	import Pkg
	Pkg.activate(Base.current_project())

	using PlutoUI
	using AlgebraOfGraphics
	using CairoMakie
	using CSV, DataFrames
	using Dates, DelimitedFiles
	using Glob
	using ImageFiltering
	using Latexify
	using LombScargle, Measurements, Statistics
end

# ╔═╡ 68f870ba-3825-4f19-a237-7f04d06ee0eb
using CommonMark

# ╔═╡ b85cbebb-3334-4672-bf36-1070fd5dff46
begin
	ENV["PYTHON"] = ""
	Pkg.build("PyCall")
	using PyCall
	const Conda = PyCall.Conda
	Conda.add("lightkurve", :WASP50b)
end

# ╔═╡ 670b88e4-1e96-494d-bfcc-235092bb6e96
md"""
# Photometric Monitoring

In this notebook we gather and analyze the available photometric data for this target.

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ 8d42d1c7-c517-41c4-9a5d-2908d2ac2463
const FIG_PATH = "figures/stellar_activity"

# ╔═╡ 0cbe4263-799f-4ee3-9a94-3ba879528b01
md"""
## $(@bind plot_ASASSN CheckBox()) ASAS-SN 🌍

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

# ╔═╡ 9094c4a4-3b75-4e21-97a7-600de734867b
describe(df_ASASSN, :all)

# ╔═╡ 2b2dd83b-ce99-4551-b8aa-79ca6db5dd06
function name(dt_julian)
	dt = julian2datetime(dt_julian)
	return "$(year(dt)) $(monthname(dt)) $(day(dt))"
end

# ╔═╡ 6eaf882c-0cb5-415f-b8fe-c071ee25a895
md"""
Looks good. According to the table above, the data spans
**from $(join(name.(extrema(df_ASASSN.hjd)), " to "))** from two cameras (**bd**, commissioned 2013 December at Haleakala; **bh**, commissioned 2015 July at CTIO), both in the V band:
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
		
	fig
end

# ╔═╡ 50ee0fbb-30f8-4e29-9de8-0173efcee364
1e6 * median(df_ASASSN.flux_err)

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
begin
	py"""
	from astropy.coordinates import SkyCoord, EarthLocation
	from astropy import units as u
	from astropy.time import Time
	
	
	# https://gist.github.com/StuartLittlefair/4ab7bb8cf21862e250be8cb25f72bb7a
	def helio_to_bary(coords, hjd, obs_name):
	    helio = Time(hjd, scale="utc", format="jd")
	    obs = EarthLocation.of_site(obs_name)
	    star = SkyCoord(coords, unit=(u.hour, u.deg)) 
	    ltt = helio.light_travel_time(star, "heliocentric", location=obs)
	    guess = helio - ltt
	    
	    # If we assume guess is correct - how far is heliocentric time away
	    # from true value?
	    delta = (
	    guess + guess.light_travel_time(star, "heliocentric", obs)
	    ).jd - helio.jd
	    
	    # Apply this correction
	    guess -= delta * u.d
	    ltt = guess.light_travel_time(star, 'barycentric', obs)
	    
	    return guess.tdb + ltt
	"""
	helio_to_bary(coords, hjd, obs_name) = py"helio_to_bary"(
		coords, hjd, obs_name
	).value
end

# ╔═╡ c86f5adc-6e17-44c4-b754-1b5c42557809
const coords_W50 = "1:12:43.2 +31:12:43" # RA, Dec of WASP-50

# ╔═╡ 5f1ca89d-61e8-4ade-93e1-13716dc5fd46
loc_to_BJD(t, coords, camera) = camera == "Haleakala" ?
	helio_to_bary(coords, t, "Haleakala") :
	helio_to_bary(coords, t, "CTIO")

# ╔═╡ 42b72605-f89e-42c4-a159-d43e4620140f
df_ASASSN.bjd = loc_to_BJD.(df_ASASSN.hjd, coords_W50, df_ASASSN.camera);

# ╔═╡ 79d08932-66f3-4ed9-bc13-f1ac3229e95d
df_ASASSN

# ╔═╡ 92548af3-9a26-4202-88f2-ba3a31181686
begin
	df_sorted = sort(df_ASASSN, :bjd)
	t_ASASSN, f_ASASSN, f_err_ASASSN = eachcol(
		df_sorted[!, [:bjd, :flux_mJy_, :flux_err]]
	)
	#t_ASASSN .-= 2.457e6
end;

# ╔═╡ 3033c5f2-dd7a-4490-9d67-0ee26d8b57a0
# lc_ASASSN = lk.LightCurve(
# 	time=t_ASASSN, flux=f_ASASSN, flux_err=f_err_ASASSN
# ).normalize()

# ╔═╡ dbe317fe-540d-44e8-b8e7-6c465c79559f
md"""
``Δt_\text{ASAS-SN}`` = $(@bind binsize_ASASSN PlutoUI.Slider(1:30, default=7, show_value=true)) days
"""

# ╔═╡ 9a195365-e68f-43e2-8870-c09153e2ce91
md"""
### Plot

With the BJD times computed, we can now plot the ASAS-SN photometry, binned to **$(binsize_ASASSN) days**:
"""

# ╔═╡ 682499cb-af79-48f7-9e74-0185894f65fe
md"""
For comparison, the average cadence for this data is $(t_ASASSN |> diff |> mean) days.
"""

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

# ╔═╡ 7370a1d9-4f8e-4788-adac-b8be2bcc9643
function plot_phot!(ax, t, f, f_err; t_offset=0.0, relative_flux=false, binsize=1.0)
	t_rel = t .- t_offset
	if relative_flux
		f_med = median(f)
		Δf, Δf_err = @. 1.0 + (f - f_med) / f_med , @. f_err / f_med
		
	else
		f_med = 1.0
		Δf, Δf_err = f, f_err
	end
	
	# Original data
	errorbars!(ax, t_rel, Δf, Δf_err, color=(:darkgrey, 0.25))
	scatter!(ax, t_rel, Δf, color=(:darkgrey, 0.25))
	
	# Binned data
	lc = lk.LightCurve(
		time=t, flux=f, flux_err=f_err
	).normalize()
	lc_binned = lc.bin(binsize).remove_nans()
	t_binned, f_binned, f_err_binned = (
		lc_binned.time.value,
		lc_binned.flux,
		lc_binned.flux_err
	)
	#bin_lc(t_rel, Δf, Δf_err, binsize)
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

# ╔═╡ 2952e971-bce5-4a1e-98eb-cb2d45c8c5a8
Time = pyimport("astropy.time").Time

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

# ╔═╡ 43de00bf-e616-43c5-92ce-1044cbd8cfe5
1e6 .* [median(lc.flux_err) for lc ∈ lcs_oot]

# ╔═╡ 2aebd994-1490-43d0-a9d1-e68652aebc32
cm"""
~~h~~
"""

# ╔═╡ 99fcf959-665b-44cf-9b5f-fd68a919f846
md"""
### $(@bind plot_pg CheckBox()) Periodogram
"""

# ╔═╡ de104bdf-e95a-4a6b-9178-c6b2a2f2f5ea
function compute_pgram(lc; min_period=0.5, max_period=30.0)
	plan = LombScargle.plan(
		lc.time.value, lc.flux .± lc.flux_err,
		minimum_frequency = 1.0 / max_period,
		maximum_frequency = 1.0 / min_period,
	)
	return lombscargle(plan), plan
end

# ╔═╡ f3425d9c-861e-4b26-b352-bd0669c7f1f9
if plot_ASASSN let
	fig = Figure(resolution=(800, 800))
	
	### Photometry plot ####
	ax = Axis(fig[1, 1], xlabel="Time (BTJD)", ylabel="Relative flux (ppm)")
	
	# Mark transit epochs
	Δjulian_transit_dates = julian_transit_dates #.- 2.457e6
	vlines!(ax, Δjulian_transit_dates;
		linestyle = :dash,
		color = :darkgrey,
	)
	
	# Label transit epochs
	# for (i, (utc, jd)) in enumerate(zip(utc_transit_dates, Δjulian_transit_dates))
	# 	text!(ax, "Transit $i\n$utc";
	# 		position = Point2f0(jd, 5.0e4),
	# 		textsize = 14,
	# 		align = (:left, :center),
	# 		offset = Point2f0(10, 0),
	# 		color = :grey,
	# 	)
	# end
	
	ax_phot, t_binned, f_binned, f_err_binned, lc_binned, f_med = plot_phot!(
		ax, t_ASASSN, f_ASASSN, f_err_ASASSN;
		t_offset=0, relative_flux=true, binsize=binsize_ASASSN
	)
	
	#### Periodogram #######
	ax_pg = Axis(fig[2, 1], xlabel="Periods (days)", ylabel="log10 Power")
	pgram, plan = compute_pgram(lc_binned)
	b = LombScargle.bootstrap(100, plan)
	P_max = findmaxperiod(pgram)[1]
	lines!(ax_pg, periodpower(pgram)...;
			label = "P_max: $(P_max) days",
		)
	hlines!(ax_pg, collect(fapinv.(Ref(b), (0.01, 0.05, 0.1))), linestyle=:dash)
	axislegend(ax_pg, position=:lc)
			
	fig
	end
end

# ╔═╡ 2215ed86-fa78-4811-88ab-e3521e4a1dea
function compute_window_func(lc; min_period=0.5, max_period=30.0)
	t = lc.time.value
	f = oneunit.(t)
	f_err = median(lc.flux_err) .* f
	lc_window_func = lk.LightCurve(time=t, flux=f, flux_err=f_err)
	return compute_pgram(lc_window_func; min_period=min_period, max_period=max_period)
end

# ╔═╡ d7f034c5-5925-4b91-9bea-1068a7ce9252
begin
	pgrams, pgrams_window, plans, P_maxs = [], [], [], []
	for lc in lcs_oot
		pgram, plan = compute_pgram(lc)
		pgram_window, _ = compute_window_func(lc)
		P_max = findmaxperiod(pgram)[1]
		
		push!(pgrams, pgram)
		push!(pgrams_window, pgram_window )
		push!(plans, plan)
		push!(P_maxs, P_max)
	end
end

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

# ╔═╡ a50ef756-ade6-48a3-8d3a-17b56ce03c26
md"""
### $(@bind plot_folded CheckBox()) Folded lightcurves
"""

# ╔═╡ 3128e57f-df4f-4811-b867-8a293d7d536d
function compute_pgram_model(lc, P)
	t_fit = lc.time.value
	s_fit = LombScargle.model(
		lc.time.value,
		lc.flux,
		lc.flux_err,
		inv(P),
	)
	lc_fit = lk.LightCurve(time=t_fit, flux=s_fit)
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
		#Δt = (lc_folded.time.value |> diff |> median) * 5
		push!(lcs_folded, lc_folded)
		push!(lcs_folded_binned, lc_folded.bin(bins=200))
		
		# Model
		lc_fit_folded = compute_pgram_model(lc, P)	
		push!(lcs_fit_folded, lc_fit_folded)
	end
end

# ╔═╡ 056281a2-4786-45eb-a9fa-57515153f66c
md"""
### Spot parameter estimation
"""

# ╔═╡ 3a612743-7071-4d85-a48d-0a4b12facffc
ΔLs = [
	minimum(lc.flux) / median(lc.flux)
	for lc in lcs_fit_folded
]

# ╔═╡ 278ea804-e2dd-4ca8-9a20-0e9a25746e02
T₀ = 5_520

# ╔═╡ d26122f1-1602-440c-8ba9-72469a782104
f_sp(T_sp, ΔL, T₀=T₀) = (T₀^4 /  (T₀^4 - T_sp^4)) * (1 - ΔL)

# ╔═╡ c0d53fc0-467b-4e49-a829-f149e14e3d08
[f_sp.((2_200, 2_800), ΔL, T₀) for ΔL in ΔLs]

# ╔═╡ 3077c3bf-9ddd-46db-94b7-b7a8120f1485
Ts = 10.0:10.0:5_000

# ╔═╡ df370404-2f12-4925-8827-6198793ae842
extrema(f_sp.(Ts, ΔLs[end], T₀))

# ╔═╡ 8bd502e8-e67d-44be-8a60-d7ad2c147d70
lcs_oot_comb = lcs_oot[end]

# ╔═╡ 3551787f-0a83-408f-9d78-41309ae3dae3
(lcs_oot_comb.flux_err |> median)

# ╔═╡ 8c7dcfab-a357-4024-94f3-42d1df80c3c2
P_maxs

# ╔═╡ 06abb8cb-9acb-49ba-81b6-37b9f52c89b1
function yee(lc)
	lcs_folded = []
	lcs_folded_binned = []
	lcs_fit_folded = []
	for P ∈ (16.3)
		# Data
		lc_folded = lc.fold(P)
		push!(lcs_folded, lc_folded)
		push!(lcs_folded_binned, lc_folded.bin(bins=200))
		
		# Model
		lc_fit_folded = compute_pgram_model(lc, P)	
		push!(lcs_fit_folded, lc_fit_folded)
	end
	
	return lcs_folded, lcs_folded_binned, lcs_fit_folded
end

# ╔═╡ 7a9dd8e0-3c2d-4c99-ae86-401554ad8558
x = yee(lcs_oot_comb)

# ╔═╡ 241d1675-007d-4b85-a818-d792681cb45a
x[3]

# ╔═╡ 2429035b-5b8e-45d5-9957-99ad772324af
ΔLs2 = [
	minimum(lc.flux) / median(lc.flux)
	for lc in x[3]
]

# ╔═╡ 377c1376-b81f-40e5-8ab3-22cc7d77d7a5
extrema(f_sp.(Ts, ΔLs2[end], T₀))

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

# ╔═╡ ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
md"""
## Notebook setup
"""

# ╔═╡ 79acbb60-803a-4047-b26d-1cf6262274a0
function savefig(fig, fpath)
	mkpath(dirname(fpath))
    save(fpath, fig)
	@info "Saved to: $(fpath)"
end

# ╔═╡ 9e2ce576-c9bd-11eb-0699-47af13e79589
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

# ╔═╡ 82222ee8-f759-499d-a072-c219cc33ccad
let
	fig = Figure(resolution=FIG_WIDE)
	
	for (i, (lc, lc_oot)) ∈ enumerate(zip(lcs_cleaned[1:end-1], lcs_oot[1:end-1]))
		ax = Axis(fig[i, 1])
		errorbars!(ax, lc.time.value, lc.flux, lc.flux_err;
			color = (:darkgrey, 0.25),
			markersize = 15,
			label = """
			Sector $(lc.meta["SECTOR"]), $(lc.meta["AUTHOR"])
			"""
		)
		
		ylims!(ax, 0.97, 1.02)
		#scatter!(fig[i, 1], lc.time.value, lc.flux)
		
		scatter!(ax, lc_oot.time.value, lc_oot.flux;
			color = :darkgrey, label="OOT baseline",
			#markersize = 5,
		)

		axislegend()
	end

	linkyaxes!(filter(x -> x isa Axis, fig.content)...)

	Label(fig[end+1, 1], "Time (BTJD days)", tellwidth=false)
	Label(fig[1:end-1, 0], "Relative flux", rotation=π/2)

	savefig(fig, "$(FIG_PATH)/TESS_flux.png")
	
	fig
end

# ╔═╡ 94d05a5b-b05e-4407-bcd3-7d625680a262
if plot_pg let
	fig = Figure(resolution=FIG_WIDE)
	
	ax_window = Axis(fig[1, 1])
	for (i, pgram_window) ∈ enumerate(pgrams_window)
		lines!(ax_window, periodpower(pgram_window)..., color=COLORS_SERIES[i])
	end
	text!(ax_window, "Window function", position=(1, 0.3))
	
	ax = Axis(fig[2, 1], xlabel="Period (days)")
	sectors = ("Sector 04", "Sector 31", "Combined")
	for (i, (pgram, plan, P_max, sector)) ∈ enumerate(
		zip(pgrams, plans, P_maxs, sectors)
	)
		# Compute FAPs
		b = LombScargle.bootstrap(100, plan)
		
		# Plot
		lines!(ax, periodpower(pgram)...;
			color = COLORS_SERIES[i],
			label="$(sector): $(round(P_max, digits=2))"
		)
		hlines!(ax, collect(fapinv.(Ref(b), (0.01, 0.05, 0.1))), color=:darkgrey)
	end
	
	axislegend("P_max (days)", position=(0.05, 0.8))
	
	hidexdecorations!(ax_window)
	linkaxes!(ax_window, ax)
	
	#xlims!(ax, 0, 10)
	#ylims!(ax, 0, 0.3)
	
	Label(fig[1:2, 0], "Normalized power", rotation=π/2)
	
	#axislegend()
	
	savefig(fig, "$(FIG_PATH)/stellar_activity_pg.png")
	
	fig
	end
end

# ╔═╡ 49bcddbe-d413-48ae-91d8-92bcebf40518
if plot_folded let
	fig = Figure(resolution=FIG_WIDE)
	
	axs = []
	sectors = ("Sector 04", "Sector 31", "Combined")
	for (i, (lc_folded, lc_folded_binned, lc_fit_folded)) in enumerate(zip(
				lcs_folded, lcs_folded_binned, lcs_fit_folded
		))
		ax = Axis(fig[i, 1])
		push!(axs, ax)
		scatter!(ax, lc_folded.time.value, lc_folded.flux, color=(:darkgrey, 0.5))
		scatter!(ax, lc_folded_binned.time.value, lc_folded_binned.flux, color=COLORS[2])
		lines!(ax, lc_fit_folded.time.value, lc_fit_folded.flux, color=0.5 .*(COLORS[2], 1.0))
		text!(ax, "$(sectors[i])";
			position = (3.8, 1.006),
		)
		ylims!(ax, 0.991, 1.01)
	end
	
	linkaxes!(axs...)
	
	axs[end].xlabel = "Phase"
	axs[2].ylabel = "Normalized flux"
	
	savefig(fig, "$(FIG_PATH)/stellar_activity_phase.png")
	
	fig
	end
end

# ╔═╡ 61c6dd34-5712-4077-abae-2ee2634dc709
@with_terminal Conda.list(:WASP50b)

# ╔═╡ Cell order:
# ╟─670b88e4-1e96-494d-bfcc-235092bb6e96
# ╟─8d42d1c7-c517-41c4-9a5d-2908d2ac2463
# ╟─0cbe4263-799f-4ee3-9a94-3ba879528b01
# ╟─b00c28a2-26b1-442e-a347-39fb66b825a0
# ╠═fa233a2c-7e89-4e71-85b8-824c5c341650
# ╠═9094c4a4-3b75-4e21-97a7-600de734867b
# ╟─6eaf882c-0cb5-415f-b8fe-c071ee25a895
# ╟─2b2dd83b-ce99-4551-b8aa-79ca6db5dd06
# ╠═98704543-8cb7-4fca-b601-2a8d2dfa4833
# ╠═8e6008ca-a762-4450-a09e-bd0c2bbac4f2
# ╠═50ee0fbb-30f8-4e29-9de8-0173efcee364
# ╟─7e5ea835-8eb2-44d5-905d-433767b6b91a
# ╟─5bc41820-3782-4f82-80b6-6da62805ca8f
# ╠═5bf1136b-2d13-4463-8d74-9ade1e2cee96
# ╠═c86f5adc-6e17-44c4-b754-1b5c42557809
# ╠═5f1ca89d-61e8-4ade-93e1-13716dc5fd46
# ╠═42b72605-f89e-42c4-a159-d43e4620140f
# ╠═79d08932-66f3-4ed9-bc13-f1ac3229e95d
# ╟─9a195365-e68f-43e2-8870-c09153e2ce91
# ╠═92548af3-9a26-4202-88f2-ba3a31181686
# ╠═3033c5f2-dd7a-4490-9d67-0ee26d8b57a0
# ╠═dbe317fe-540d-44e8-b8e7-6c465c79559f
# ╠═f3425d9c-861e-4b26-b352-bd0669c7f1f9
# ╠═7370a1d9-4f8e-4788-adac-b8be2bcc9643
# ╟─682499cb-af79-48f7-9e74-0185894f65fe
# ╟─78d85c7c-da06-41ab-915b-48d93a010967
# ╟─97e7feee-11b2-4a35-9327-b5c0d05b2a23
# ╠═0c790d2f-64d4-4e13-9629-a9725cd7086d
# ╠═2952e971-bce5-4a1e-98eb-cb2d45c8c5a8
# ╠═ec12acb8-9124-4cc0-8c9f-6525c1565dfd
# ╠═708d54a5-95fd-4f15-9681-f6d8e7b9b05c
# ╟─34fcd73d-a49c-4597-8e63-cfe2495eee48
# ╠═dff46359-7aec-4fa1-bc7a-89785dfca0e8
# ╠═6e62b3cc-96ce-43fd-811b-4b2b102cfd61
# ╟─241c462c-3cd9-402d-b948-b9b1f608b727
# ╠═31d5bc92-a1f2-4c82-82f2-67755f9aa235
# ╠═82222ee8-f759-499d-a072-c219cc33ccad
# ╠═3551787f-0a83-408f-9d78-41309ae3dae3
# ╠═43de00bf-e616-43c5-92ce-1044cbd8cfe5
# ╠═68f870ba-3825-4f19-a237-7f04d06ee0eb
# ╠═2aebd994-1490-43d0-a9d1-e68652aebc32
# ╟─99fcf959-665b-44cf-9b5f-fd68a919f846
# ╠═94d05a5b-b05e-4407-bcd3-7d625680a262
# ╠═d7f034c5-5925-4b91-9bea-1068a7ce9252
# ╠═de104bdf-e95a-4a6b-9178-c6b2a2f2f5ea
# ╠═2215ed86-fa78-4811-88ab-e3521e4a1dea
# ╠═d1f7ed4b-4599-48bd-aac5-93920dae9151
# ╟─a50ef756-ade6-48a3-8d3a-17b56ce03c26
# ╠═49bcddbe-d413-48ae-91d8-92bcebf40518
# ╠═97ced6ba-ff74-46b4-90d5-18e7b2f1b903
# ╠═3128e57f-df4f-4811-b867-8a293d7d536d
# ╟─056281a2-4786-45eb-a9fa-57515153f66c
# ╠═3a612743-7071-4d85-a48d-0a4b12facffc
# ╠═278ea804-e2dd-4ca8-9a20-0e9a25746e02
# ╠═d26122f1-1602-440c-8ba9-72469a782104
# ╠═c0d53fc0-467b-4e49-a829-f149e14e3d08
# ╠═3077c3bf-9ddd-46db-94b7-b7a8120f1485
# ╠═df370404-2f12-4925-8827-6198793ae842
# ╠═377c1376-b81f-40e5-8ab3-22cc7d77d7a5
# ╠═8bd502e8-e67d-44be-8a60-d7ad2c147d70
# ╠═8c7dcfab-a357-4024-94f3-42d1df80c3c2
# ╠═06abb8cb-9acb-49ba-81b6-37b9f52c89b1
# ╠═7a9dd8e0-3c2d-4c99-ae86-401554ad8558
# ╠═241d1675-007d-4b85-a818-d792681cb45a
# ╠═2429035b-5b8e-45d5-9957-99ad772324af
# ╟─18223d42-66d8-40d1-9d89-be8af46853e2
# ╠═682c3732-e68f-4fdb-bd63-553223308364
# ╟─ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
# ╟─79acbb60-803a-4047-b26d-1cf6262274a0
# ╠═9e2ce576-c9bd-11eb-0699-47af13e79589
# ╠═61c6dd34-5712-4077-abae-2ee2634dc709
# ╠═b85cbebb-3334-4672-bf36-1070fd5dff46
# ╠═55beac98-0929-4a55-91f7-cee7c781498c
