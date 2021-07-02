import Pluto

names = [
    "raw_frames",
    "reduced",
    "reduced_LDSS3",
    "detrended_wlcs",
    "detrended_bwlcs",
    "transmission_spectra",
    "photometric_monitoring",
    "photometric_monitoring2",
]

s = Pluto.ServerSession();
for name in names
    println("\nRunning $(name).jl ...")
    nb = Pluto.SessionActions.open(s, "notebooks/$(name).jl"; run_async=false)
    html_contents = Pluto.generate_html(nb)
    write("html/$(name).html", html_contents)
end
