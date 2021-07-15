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
	using Measurements
	using Measurements: value, uncertainty
	using NaturalSort
	using OrderedCollections
	using Printf
	using Statistics
	using PlutoUI: TableOfContents, Select, Slider, as_svg, with_terminal
	using Unitful
	
	# Python setup
	ENV["PYTHON"] = "/home/mango/miniconda3/envs/WASP50/bin/python"
	Pkg.build("PyCall")
	using PyCall
	
	##############
	# PLOT CONFIGS
	##############
	const FIG_TALL = (900, 1_200)
	const FIG_WIDE = (1_350, 800)
	
	# set_aog_theme!()
	# update_theme!(
	# 	Theme(
	# 		Axis = (xlabelsize=18, ylabelsize=18,),
	# 		Label = (textsize=18,  padding=(0, 10, 0, 0)),
	# 		Lines = (linewidth=3, cycle=Cycle([:color, :linestyle], covary=true)),
	# 		Scatter = (linewidth=10,),
	# 		fontsize = 18,
	# 		rowgap = 0,
	# 		colgap = 0,
	# 	)
	# )
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
df_ASASSN[!, :bjd] = loc_to_BJD.(df_ASASSN.hjd, coords_W50, df_ASASSN.camera);

# ╔═╡ 79d08932-66f3-4ed9-bc13-f1ac3229e95d
df_ASASSN

# ╔═╡ ebf9fef2-8be0-4993-92c9-4b661741c380
let
	t = df_ASASSN.bjd .- 2.457e6
	f = df_ASASSN.flux_mJy_
	f_norm = (f .- mean(f)) ./ f .+ 1.0
	f_err = df_ASASSN.flux_err ./ f
	#errorbars(t, f_norm, f_err)
	
	CSV.write(
		"/home/mango/Desktop/WASP50LC_ASASSN.txt",
		DataFrame(:t=>t, :f =>f_norm, :f_err=>f_err),
		delim = "\t",
	)		
end

# ╔═╡ 92548af3-9a26-4202-88f2-ba3a31181686
begin
	df_sorted = sort(df_ASASSN, :bjd)
	t_ASASSN, f_ASASSN, f_err_ASASSN = eachcol(df_sorted[!, [:bjd, :mag, :mag_err]])
end;

# ╔═╡ dbe317fe-540d-44e8-b8e7-6c465c79559f
md"""
``t_\text{window}`` = $(@bind t_window Slider(1:20, default=7, show_value=true)) days
"""

# ╔═╡ 9a195365-e68f-43e2-8870-c09153e2ce91
md"""
### Plot

With the BJD times computed, we can now plot the ASAS-SN photometry, binned to **$(t_window) days**:
"""

# ╔═╡ da96bd4e-b755-4577-9f2f-783a3f5d619f
md"""
!!! warning
	Put back BTJD
"""

# ╔═╡ 682499cb-af79-48f7-9e74-0185894f65fe
md"""
For comparison, the cadence for this data has an average of $(t_ASASSN |> diff |> mean) days.
"""

# ╔═╡ 78d85c7c-da06-41ab-915b-48d93a010967
md"""
## TESS 🌌
We next turn to the TESS photometry, before finally combining with the ASAS-SN data.
"""

# ╔═╡ 97e7feee-11b2-4a35-9327-b5c0d05b2a23
md"""
### Light curve collection

First we use [`lightkurve`](https://docs.lightkurve.org/whats-new-v2.html) to download the available data products from TESS:
"""

# ╔═╡ 7952b639-d10f-4d2d-9289-4b92b2d1c60d
begin
	py"""
	import transitleastsquares as tls
	import lightkurve as lk
	
	def search_results(target, kwargs):
		return lk.search_lightcurve(target, **kwargs)
	"""
	search_results(target; kwargs...) = py"search_results"(target; kwargs)
end

# ╔═╡ cc76af35-a066-4d42-be49-eefec91f7d96
all_srs = search_results("WASP-50")

# ╔═╡ dff46359-7aec-4fa1-bc7a-89785dfca0e8
srs = search_results("WASP-50"; author=["TESS-SPOC"])

# ╔═╡ 34fcd73d-a49c-4597-8e63-cfe2495eee48
md"""
From the $(length(all_srs)) data products found, we see that $(length(srs)) are available from the [Science Processing Operations Center (SPOC)](https://heasarc.gsfc.nasa.gov/docs/tess/pipeline.html):
"""

# ╔═╡ 6e62b3cc-96ce-43fd-811b-4b2b102cfd61
lcs = srs.download_all();

# ╔═╡ 241c462c-3cd9-402d-b948-b9b1f608b727
md"""
We show the normalized PDCSAP flux below for each sector: 
"""

# ╔═╡ 82222ee8-f759-499d-a072-c219cc33ccad
let
	fig = Figure(resolution=FIG_TALL, title="yee")
		
	for (i, lc) in enumerate(lcs)
		lines(
			fig[i, 1],
			lc.normalize().remove_nans().time.value,
			lc.normalize().remove_nans().flux,
			label = """
			Sector $(lc.meta["SECTOR"]), $(lc.meta["AUTHOR"])
			$(lc.meta["EXPOSURE"]) s
			$(srs[i].exptime)
			"""
		)
		# lines!(
		# 	fig[i, 1],
		# 	lc.normalize().remove_nans().time.value,
		# 	lc.normalize().remove_nans().flux.value,
		# )

		axislegend()
	end

	#linkyaxes!(filter(x -> x isa Axis, fig.content)...)

	Label(fig[end+1, 1], "Time (BTJD days)", tellwidth=false)
	Label(fig[1:end-1, 0], "Flux (e⁻/s)", rotation=π/2)
	
	fig
end

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

# ╔═╡ 3327596c-56f1-4024-9490-ee69bd514007
md"""
### Light curve binning
"""

# ╔═╡ e7eba854-3846-4ba4-bbdb-c5f7dbdf08a7
df_TESS_S04 = CSV.File(
	"data/photometric/TESS_baseline_sector_04.csv",
) |> DataFrame

# ╔═╡ e1648b37-b52b-4c46-a5fb-4e5bc2743cd4
df_TESS_S31 = CSV.File(
	"data/photometric/TESS_baseline_sector_31.csv",
) |> DataFrame

# ╔═╡ 09c666f7-a8b4-4b47-b9fd-1351e8bd28f9
t_window_TESS = 200.0 * (1/86_400) * 10

# ╔═╡ 5cca08d9-e035-4a15-b575-069e6b89b6db
t_TESS_S04, f_TESS_S04, f_err_TESS_S04 = eachcol(
	df_TESS_S04[!, [:time, :flux, :flux_err]]
)

# ╔═╡ 8da5f6ed-f0d2-484e-a5ae-4a6ac603f0df
(t_TESS_S04[end] + 2.457e6) |> julian2datetime

# ╔═╡ a87d6163-4fd9-49c3-a83b-ab9c727edf99
t_TESS_S31, f_TESS_S31, f_err_TESS_S31 = eachcol(
	df_TESS_S31[!, [:time, :flux, :flux_err]]
)

# ╔═╡ fe062294-c575-451c-935d-0d75194b5137
const TESSMAG = 11.05290031

# ╔═╡ faff66de-c844-49be-81a0-d019780a49ed
flux_to_mag(f, M=TESSMAG) = M - 2.5*log10(f)

# ╔═╡ ea8d87e4-4254-4e4e-a6ab-a2ef39d9a756
flux_to_mag_err(f) = -2.5*log10(1.0 - f)

# ╔═╡ 99fcf959-665b-44cf-9b5f-fd68a919f846
md"""
### Periodogram
"""

# ╔═╡ ec12acb8-9124-4cc0-8c9f-6525c1565dfd
begin
	py"""
	import lightkurve as lk
	import astropy.units as u
	from tess_sip import SIP
	
	def run_sip(target, srs_kwargs, sip_kwargs):
		srs = lk.search_lightcurve(target, **srs_kwargs)
		lcs = srs.download_all()
		r = SIP(lcs, **sip_kwargs)
		return r, srs
	"""
	run_sip(target; srs_kwargs, sip_kwargs) = py"run_sip"(
		target; srs_kwargs, sip_kwargs
	)
end

# ╔═╡ 8f382cc4-2c3b-4d40-83e7-044fbb6efb2e
sip_kwargs = Dict(
	"max_period" => 30.0,
	"nperiods" => 400,
	"bin_kwargs" => Dict("time_bin_size"=>(12.0u"hr" |> u"d").val),
)

# ╔═╡ f287c0a6-e5d9-43be-b0b9-ded9273bdfc1
r_S04, srs_S04 = run_sip(
	"WASP50",
	srs_kwargs = Dict("sector"=>4, "author"=>"TESS-SPOC"),
	sip_kwargs = sip_kwargs
); srs_S04

# ╔═╡ 4c078bc4-b05d-4357-8dc1-d660b35cb2e0
r_S31, srs_S31 = run_sip(
	"WASP50",
	srs_kwargs = Dict("sector"=>31, "author"=>"TESS-SPOC"),
	sip_kwargs = sip_kwargs
); srs_S31

# ╔═╡ d2ed4e4b-edc1-4dfd-adeb-c6e5d5fab923
r_combined, srs_combined = run_sip(
	"WASP50",
	srs_kwargs = Dict("author"=>"TESS-SPOC"),
	sip_kwargs = sip_kwargs
); srs_combined

# ╔═╡ fd3ffd0c-4d08-4862-bc90-a6b748e1faab
data_dict = OrderedDict(
	"Sector 04" => r_S04,
	"Sector 31" => r_S31,
	"Combined" => r_combined,
);

# ╔═╡ 4fcb3ea1-bec6-4c30-82a6-4689077cdc14
function plot_periodogram!(ax, r, P)
	periods = r["periods"]
    power_source = r["power"]
    power_background = r["power_bkg"]
    power_relative = power_source ./ power_background

	lines!(ax, periods, power_source; label="source")
	lines!(ax, periods, power_background; label="background")
	lines!(ax, periods, power_relative; label="divided"
	)
	
	P_max = periods[argmax(power_relative)]
	vlines!(ax, P, linestyle=:dash, color=:darkgrey, label="P_rot (Gillon+ 2011)")
	vlines!(ax, P_max, color=:darkgrey, label="P_max")
	text!(ax, "P_max = $(round(P_max, digits=2)) days";
		position = (30, 2),
		textsize = 16,
		align = (:right, :center),
		color = :darkgrey,
	)
end

# ╔═╡ ed93fd1f-79bd-4c5a-b822-680c6f6591d6
let
	fig = Figure(resolution=FIG_TALL)
	
	P = 16.3
	
	axs = []
	for (i, (sector, r)) in enumerate(data_dict)
		ax = Axis(fig[i, 1], yscale=log10)
		plot_periodogram!(ax, r, P)
		text!(ax, sector;
			position = (30, 4.5),
			textsize = 16,
			align = (:right, :center),
			color = :darkgrey,
		)
		push!(axs, ax)
	end
	
	linkaxes!(axs...)
	hidexdecorations!.(axs[1:2])
	
	Legend(fig[1:end, 2], axs[1], margin=(10, 0, 0, 0))
	Label(fig[1:end, 0], "Power", rotation=π/2)
	
	fig
end

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

# ╔═╡ 976bd30f-2bae-4c5d-93dc-70b4df303840
md"""
We produced the above plots by performing a rolling average over sub-groups of the data. Each sub-group contains timeseries points within $Δt$ of the first point in that sub-group. We then move on to the next group and perform the same grouping operation:

!!! note
	TODO: Find a good Δt to use
"""

# ╔═╡ 7370a1d9-4f8e-4788-adac-b8be2bcc9643
function plot_phot!(ax, t_window, t, f, f_err; t_offset=0.0, relative_flux=false)
	t_rel = t .- t_offset
	
	if relative_flux
		Δf, Δf_err = f .- mean(f), f_err
		
	else
		Δf, Δf_err = f, f_err
	end
	
	# Original data
	errorbars!(ax, t_rel, Δf, Δf_err, color=(:darkgrey, 0.25))
	scatter!(ax, t_rel, Δf, color=(:darkgrey, 0.25))
	
	# Averaged data
	t_avg, f_avg = bin_data(t_rel, Δf .± Δf_err, Δ=t_window)
	f_unc_avg = uncertainty.(f_avg) |> mean
	errorbars!(ax, t_avg, f_avg)
	scatter!(ax, t_avg, f_avg, label="avg err: $(f_unc_avg)")
	
	axislegend(ax)
	
	return ax, t_avg, value.(f_avg), uncertainty.(f_avg)
end

# ╔═╡ 2ee9b0ca-f6a2-47db-abfc-9a44e64b2e42
let
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="Time (BTJD)", ylabel="ΔMag")
	
	# Mark transit epochs
	Δjulian_transit_dates = julian_transit_dates .- 2.45e6
	vlines!(ax, Δjulian_transit_dates;
		linestyle = :dash,
		color = :darkgrey,
	)
	
	# Label transit epochs
	for (i, (utc, jd)) in enumerate(zip(utc_transit_dates, Δjulian_transit_dates))
		text!(ax, "Transit $i\n$utc";
			position = Point2f0(jd, 0.08),
			textsize = 14,
			align = (:left, :center),
			offset = Point2f0(10, 0),
			color = :grey,
		)
	end
	
	# Plot photometry measurements
	ax, t, f, f_err = plot_phot!(
		ax, t_window, t_ASASSN, f_ASASSN, f_err_ASASSN;
		t_offset=2.45e6, relative_flux=true,
	)
	
	CSV.write(
		"/home/mango/Desktop/WASP50LC_ASASSN_binned.txt",
		DataFrame(:t=>t, :f =>f, :f_err=>f_err),
		delim = "\t",
	)
end

# ╔═╡ 8906a2a2-65c9-4dc1-aaef-078a6ddaaff2
let
	fig = Figure()
	
	ax_top = Axis(fig[1, 1],)
	ax_bottom = Axis(fig[2, 1])

	a_ax, t, f, f_err = plot_phot!(
		ax_top, t_window_TESS, t_TESS_S04,
		f_TESS_S04, f_err_TESS_S04;
		relative_flux=true,
	)
	CSV.write(
		"/home/mango/Desktop/WASP50LC_S04_binned.txt",
		DataFrame(:t=>t, :f=>f, :f_err=>f_err),
		delim = "\t",
	)
	
	b_ax, t, f, f_err = plot_phot!(
		ax_bottom, t_window_TESS, t_TESS_S31,
		f_TESS_S31, f_err_TESS_S31;
		relative_flux=false,
	)
	CSV.write(
		"/home/mango/Desktop/WASP50LC_S31_binned.txt",
		DataFrame(:t=>t, :f=>f, :f_err=>f_err),
		delim = "\t",
	)
	
	fig #|> as_svg
end

# ╔═╡ ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
md"""
## Notebook setup
"""

# ╔═╡ f8e37bb8-bdd8-4eea-82c3-1b1f3a5617a1
html"""
<style>
#launch_binder {
	display: none;
}
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

# ╔═╡ 23e9940b-15fe-4120-9a17-1b3dcd043441
md"""
!!! note
	Can remove this after `tess_sip` is updated
"""

# ╔═╡ 210abd28-0759-47c2-88ed-b5bf86a004e5
begin
	py"""
	import sys
	import os
	
	def update_sys():
		sys.path.insert(0, os.getcwd())
	"""
	update_sys() = py"update_sys"()
end

# ╔═╡ 8aef8c07-51db-4004-b9af-64484f488b0b
function corner_plot!(fig; n_params=3)
	for j in 1:n_params, i in 1:n_params
		# Create subplot apply global settings
		ax = Axis(fig[i, j];
			aspect = 1,
			# xticklabelrotation = π/4,
			# xticks = LinearTicks(3),
			# yticks = LinearTicks(3),
			# limits = ((-14, 14), (-14, 14)),
		)
		# Hide upper triangle
		j > i && (hidedecorations!(ax); hidespines!(ax))
		# Hide y ticks on diagonals
		j == i && (hideydecorations!(ax); ylims!(0, 0.4))
		# Hide x ticks on all diagonal elements except the bottom one
		j == i && i != n_params && (
			hidexdecorations!(ax, grid=false);
		)
		# Hide ticks on interior lower triangle
		j < i && i != n_params && j != 1 && (
			hideydecorations!(ax, grid=false);
			hidexdecorations!(ax, grid=false);
		)
		# Hide remaining xyticks
		j < i && j == 1 && i != n_params && (
			hidexdecorations!(ax, grid=false);
		)
		j < i && i == n_params && j != 1 && (
			hideydecorations!(ax, grid=false);
		)
	end
			
	# Plot corners from each night
	#for (i, (transit, cube)) in enumerate(cubes)		
	#end
	
	# Align axes limits and apply labels
	# axs = reshape(copy(fig.content), n_params, n_params)
	# [linkxaxes!(reverse(axs[:, j])...) for j in 1:n_params]
	# for (j, (param, param_latex)) in enumerate(PARAMS)
	# 	axs[end, j].xlabel = "Δ"*param_latex
	# 	axs[j, begin].ylabel = "Δ"*param_latex
	# end
	
	# Legend(fig[1, end], elems, elem_labels;
	# 	halign = :right,
	# 	valign = :top,
	# 	patchsize = (25, 25),
	# 	tellwidth = false,
	# 	tellheight = false,
	# 	rowgap = 10,
	# 	labelsize = 25,
	# 	nbanks = 2,
	# 	orientation = :horizontal,
	# )
	
	#fig #|> as_svg
end

# ╔═╡ be2db165-b5b9-43f2-940b-0add949927d2
begin
	# Create empty corner plot grid
	fig = Figure(resolution=(1_400, 1_400))
	
	corner_plot!(fig; n_params=3)
	
	fig
end

# ╔═╡ a0489258-8e71-4677-84a4-d67fe6a2b9c6
samples = Dict(
	"params_$i" => rand(100)
	for i in 1:3
)

# ╔═╡ 00bd40fc-89b1-4f10-a433-fb12720950c3
function plot_corner!(fig, samples, params; color=:blue)
	for (j, p1) in enumerate(params), (i, p2) in enumerate(params)
		# 1D plot
		if i == j
			histogram!(fig[i, i], samples[p1])
			
			# density!(fig[i, i], samples[p1];
			# 	color = (color, 0.125),
			# 	strokewidth = 3,
			# 	strokecolor = color,
			# 	strokearound = true,
			# )
		end
		# 2D plot
		if i > j
			Z = kde((samples[p1], samples[p2]), npoints=(2^4, 2^4),)
			contour!(fig[i, j], Z;
				# levels = levels(Z.density, n_levels),
				# color = color,
				# linewidth = 3,
			)
		end
	end
end

# ╔═╡ 2adec440-0422-4dc9-80ef-ba7a3ebdf847


# ╔═╡ Cell order:
# ╟─670b88e4-1e96-494d-bfcc-235092bb6e96
# ╟─0cbe4263-799f-4ee3-9a94-3ba879528b01
# ╟─b00c28a2-26b1-442e-a347-39fb66b825a0
# ╠═fa233a2c-7e89-4e71-85b8-824c5c341650
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
# ╠═ebf9fef2-8be0-4993-92c9-4b661741c380
# ╟─9a195365-e68f-43e2-8870-c09153e2ce91
# ╠═92548af3-9a26-4202-88f2-ba3a31181686
# ╟─dbe317fe-540d-44e8-b8e7-6c465c79559f
# ╠═2ee9b0ca-f6a2-47db-abfc-9a44e64b2e42
# ╟─da96bd4e-b755-4577-9f2f-783a3f5d619f
# ╟─682499cb-af79-48f7-9e74-0185894f65fe
# ╟─78d85c7c-da06-41ab-915b-48d93a010967
# ╟─97e7feee-11b2-4a35-9327-b5c0d05b2a23
# ╠═7952b639-d10f-4d2d-9289-4b92b2d1c60d
# ╠═cc76af35-a066-4d42-be49-eefec91f7d96
# ╟─34fcd73d-a49c-4597-8e63-cfe2495eee48
# ╠═dff46359-7aec-4fa1-bc7a-89785dfca0e8
# ╠═6e62b3cc-96ce-43fd-811b-4b2b102cfd61
# ╟─241c462c-3cd9-402d-b948-b9b1f608b727
# ╠═82222ee8-f759-499d-a072-c219cc33ccad
# ╠═bf01b1c2-e4c6-4bed-8e79-a20f273cc387
# ╠═8da5f6ed-f0d2-484e-a5ae-4a6ac603f0df
# ╟─3327596c-56f1-4024-9490-ee69bd514007
# ╠═e7eba854-3846-4ba4-bbdb-c5f7dbdf08a7
# ╠═e1648b37-b52b-4c46-a5fb-4e5bc2743cd4
# ╠═09c666f7-a8b4-4b47-b9fd-1351e8bd28f9
# ╠═5cca08d9-e035-4a15-b575-069e6b89b6db
# ╠═a87d6163-4fd9-49c3-a83b-ab9c727edf99
# ╠═8906a2a2-65c9-4dc1-aaef-078a6ddaaff2
# ╠═fe062294-c575-451c-935d-0d75194b5137
# ╠═faff66de-c844-49be-81a0-d019780a49ed
# ╠═ea8d87e4-4254-4e4e-a6ab-a2ef39d9a756
# ╟─99fcf959-665b-44cf-9b5f-fd68a919f846
# ╠═ec12acb8-9124-4cc0-8c9f-6525c1565dfd
# ╠═f287c0a6-e5d9-43be-b0b9-ded9273bdfc1
# ╠═4c078bc4-b05d-4357-8dc1-d660b35cb2e0
# ╠═d2ed4e4b-edc1-4dfd-adeb-c6e5d5fab923
# ╠═fd3ffd0c-4d08-4862-bc90-a6b748e1faab
# ╠═8f382cc4-2c3b-4d40-83e7-044fbb6efb2e
# ╠═ed93fd1f-79bd-4c5a-b822-680c6f6591d6
# ╠═4fcb3ea1-bec6-4c30-82a6-4689077cdc14
# ╟─18223d42-66d8-40d1-9d89-be8af46853e2
# ╠═682c3732-e68f-4fdb-bd63-553223308364
# ╟─976bd30f-2bae-4c5d-93dc-70b4df303840
# ╠═7370a1d9-4f8e-4788-adac-b8be2bcc9643
# ╟─ded3b271-6b4e-4e68-b2f6-fa8cfd52c0bd
# ╠═9e2ce576-c9bd-11eb-0699-47af13e79589
# ╟─f8e37bb8-bdd8-4eea-82c3-1b1f3a5617a1
# ╟─23e9940b-15fe-4120-9a17-1b3dcd043441
# ╠═210abd28-0759-47c2-88ed-b5bf86a004e5
# ╠═be2db165-b5b9-43f2-940b-0add949927d2
# ╠═8aef8c07-51db-4004-b9af-64484f488b0b
# ╠═a0489258-8e71-4677-84a4-d67fe6a2b9c6
# ╠═00bd40fc-89b1-4f10-a433-fb12720950c3
# ╠═2adec440-0422-4dc9-80ef-ba7a3ebdf847
