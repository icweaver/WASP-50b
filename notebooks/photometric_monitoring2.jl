### A Pluto.jl notebook ###
# v0.14.8

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

# ╔═╡ fbed4840-cada-11eb-0c6c-d38d9654159f
using Plots, TimeSeries, CSV, DataFrames, PlutoUI, Statistics

# ╔═╡ 6fd6be4c-0de1-4f12-8cd3-b773eb3666f2
md"""
# Photometric Monitoring ⭐

$(TableOfContents(title="📖 Table of Contents"))
"""

# ╔═╡ d8f23d33-3ac3-417b-a2f5-cf7dc61beb01
md"""
## Load data
"""

# ╔═╡ 3cec41d5-bb87-446b-b39e-80b949d108da
md"""
Let's start by loading in the data and summarizing its contents:
"""

# ╔═╡ 4834702f-7dd6-4e9b-bfb7-8c9c2ba7651b
df_phot_mon = CSV.File(
	"data/photometric/AP37847073.csv",
	normalizenames = true,
) |> DataFrame

# ╔═╡ ed9ea9e6-9728-4313-b082-7d8515a24591
describe(df_phot_mon, :all)

# ╔═╡ 771bb36b-affa-491b-b732-2bfabd76d8ef
function name(dt_hjd)
	dt = julian2datetime(dt_hjd)
	return "$(year(dt)) $(monthname(dt)) $(day(dt))"
end

# ╔═╡ 04761ee0-2f91-4b4c-9100-412f150d3cf2
md"""
According to the table above, the data spans from
$(join(name.(extrema(df_phot_mon.hjd)), " - ")) from two cameras (bd: 2013 December, bh: 2015 July), both in the V band. Since we will be using both cameras, let's sort the entire dataframe by time and plot the results next.
"""

# ╔═╡ 533d8696-425a-4788-95ef-1370d548ada3
md"""
## Plot
"""

# ╔═╡ 1d789b40-8fa1-4983-aac2-7275e5a0e29c
md"""
We show the original data in blue, and rolling averaged data in orange
"""

# ╔═╡ 8fca400e-09bc-4fb8-8219-3450bf1e16d4
plotly() # Use the interactive plotly backend

# ╔═╡ 7cb85f49-b3c6-4145-add2-6719e2c7201e
@bind window Slider(1:10, show_value=true)

# ╔═╡ 126a6afa-6959-4a88-bb19-d37249468778
begin
	# Sort data and convert HJD -> UTC
	df_sorted = sort(df_phot_mon, :hjd)
	t, f, f_err = eachcol(df_sorted[!, [:hjd, :mag, :mag_err]])
	t = julian2datetime.(t)
	
	# Compute rolling average
	ta = moving(median, TimeArray(t, f), window)

	# Plot
	p = scatter(t, f; yerr=f_err, xlabel="Time (HJD)", ylabel="Magnitude")
	scatter!(p, ta)
end

# ╔═╡ cba385b1-7d69-408a-a622-97e6498ae978
md"""
## Packages
"""

# ╔═╡ Cell order:
# ╟─6fd6be4c-0de1-4f12-8cd3-b773eb3666f2
# ╟─d8f23d33-3ac3-417b-a2f5-cf7dc61beb01
# ╟─3cec41d5-bb87-446b-b39e-80b949d108da
# ╠═4834702f-7dd6-4e9b-bfb7-8c9c2ba7651b
# ╠═ed9ea9e6-9728-4313-b082-7d8515a24591
# ╟─04761ee0-2f91-4b4c-9100-412f150d3cf2
# ╠═771bb36b-affa-491b-b732-2bfabd76d8ef
# ╟─533d8696-425a-4788-95ef-1370d548ada3
# ╟─1d789b40-8fa1-4983-aac2-7275e5a0e29c
# ╠═8fca400e-09bc-4fb8-8219-3450bf1e16d4
# ╠═7cb85f49-b3c6-4145-add2-6719e2c7201e
# ╠═126a6afa-6959-4a88-bb19-d37249468778
# ╟─cba385b1-7d69-408a-a622-97e6498ae978
# ╠═fbed4840-cada-11eb-0c6c-d38d9654159f
