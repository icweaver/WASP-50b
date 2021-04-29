### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ 1dc0c992-a894-11eb-2eab-f1fccba239f0
using PlutoUI, CairoMakie, Colors, CSV, DataFrames, DataFramesMeta, Glob, Measurements, Statistics

# ╔═╡ c53be9cf-7722-4b43-928a-33e7b0463330
const DATA_DIR = "data/detrended/out_l/WASP50"

# ╔═╡ 5c4fcb25-9a26-43f1-838b-338b33fb9ee6
cube = Dict(
	"tspec" => Dict(
		"Transit $i" => CSV.File(fpath) |> DataFrame
		for (i, fpath) in enumerate(sort(glob("$(DATA_DIR)/w50_*/transpec.csv")))
	),

	"BMA_WLC" => Dict(
		"Transit $i" => CSV.File(
			fpath,
			comment = "#",
			delim = ' ',			
			ignorerepeated = true,
		) |> DataFrame
		for (i, fpath) in enumerate(sort(glob(
			"$(DATA_DIR)/w50_*/white-light/results.dat"
		)))
	),
)

# ╔═╡ 0ffded2f-eefe-43e0-8d56-dc96f676f32f
const N_NIGHTS = length(cube["tspec"])

# ╔═╡ 3133ffef-decc-4d6d-a317-7c7ccb2b561d
const N_BINS = length(cube["tspec"]["Transit 3"][!, "Depth (ppm)"])

# ╔═╡ 3ba614f2-601c-4381-80e0-2d714db5505e
δ(p) = p^2 * 1e6

# ╔═╡ bef0918c-c645-4557-a2e5-00b6c26573bc
const COLORS =  parse.(Colorant,
	[
		"#5daed9",  # Cyan
		"plum",
		"#f7ad4d",  # Yellow
		"mediumaquamarine",
		"#126399",  # Blue
		"#956cb4",  # Purple
		"#ff7f00",  # Orange
		"#029e73",  # Green
		"slategray",
	]
);

# ╔═╡ 855095aa-c7e1-4799-9fcb-070a95bf7656
begin
	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="hey")
	
	tspecs = Matrix{Measurement}(undef, N_BINS, N_NIGHTS) # Wavelength binned fluxes
	ps = Matrix{Measurement}(undef, 1, N_NIGHTS) # Rₚ/Rₛ
	for (i, (transit, df)) in enumerate(cube["tspec"])
		wav, depth, depth_u, depth_d = eachcol(
			df[!, ["Wav_d", "Depth (ppm)", "Depthup (ppm)", "DepthDown (ppm)"]]
		)
		max_depth_err = maximum(hcat(depth_u, depth_d), dims=2) |> vec
		
		errorbars!(ax, wav, depth, max_depth_err;
			color = COLORS[i],
			whiskerwidth = 0,
			linewidth = 3,
		)
		scatter!(ax, wav, depth;
			color = COLORS[i],
			strokewidth = 0,
			markersize = 12,
			label = transit,
		)
		
		if transit == "Transit 1"
			#push!(tspecs, depth[5:end] .± max_depth_err[5:end])
			tspecs[:, i] .= depth[5:end] .± max_depth_err[5:end]
		else
			#push!(tspecs, depth .± max_depth_err)
			tspecs[:, i] .= depth .± max_depth_err
		end
		
		symbol, p, p_u, p_d = eachcol(
			@where(cube["BMA_WLC"][transit], :Variable .== "p")
		)
		#push!(ps, p[1] ± maximum((p_u[1], p_d[1])))
		ps[1, i] = p[1] ± maximum((p_u[1], p_d[1]))
	end
	
	axislegend(orientation=:horizontal)
	
	fig
end

# ╔═╡ 571594bd-3750-4ba2-89cd-dc362f67e8b9
Measurements.value.(tspecs[:, 1])

# ╔═╡ fdf3855a-4a10-4a54-b0bb-7619111c51d6
wlc_depths = δ.(ps)

# ╔═╡ 520c1729-3c01-4cd1-adce-26e209969bca
offsets = (wlc_depths .- mean(wlc_depths))

# ╔═╡ cbc5000d-9413-4fa7-9235-c714d4682eac
tspecs

# ╔═╡ e3fc5b28-400f-43e8-a2ab-d910a030fd43
function tspec_plot(ax, X::Matrix{Measurement})
	for (i, tspec) in enumerate(eachcol(X))
		scatter!(ax, Measurements.value.(tspec), color=COLORS[i])
	end
end

# ╔═╡ 161fa635-1944-4a5f-b74a-968662c8d21d
let
	fig = Figure()
	ax = Axis(fig[1, 1])
	
	tspec_plot(ax, tspecs[:, 1:2])
	
	fig
end

# ╔═╡ Cell order:
# ╠═1dc0c992-a894-11eb-2eab-f1fccba239f0
# ╠═c53be9cf-7722-4b43-928a-33e7b0463330
# ╠═5c4fcb25-9a26-43f1-838b-338b33fb9ee6
# ╠═0ffded2f-eefe-43e0-8d56-dc96f676f32f
# ╠═3133ffef-decc-4d6d-a317-7c7ccb2b561d
# ╠═855095aa-c7e1-4799-9fcb-070a95bf7656
# ╠═e3fc5b28-400f-43e8-a2ab-d910a030fd43
# ╠═161fa635-1944-4a5f-b74a-968662c8d21d
# ╠═571594bd-3750-4ba2-89cd-dc362f67e8b9
# ╠═fdf3855a-4a10-4a54-b0bb-7619111c51d6
# ╠═520c1729-3c01-4cd1-adce-26e209969bca
# ╠═cbc5000d-9413-4fa7-9235-c714d4682eac
# ╠═3ba614f2-601c-4381-80e0-2d714db5505e
# ╠═bef0918c-c645-4557-a2e5-00b6c26573bc
