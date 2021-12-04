# This script checks the mask coordinates obtained with get_mask_coords.py.
#
# This chips are plotted as:
# 1 2 3 4
# 6 5 8 7
#
# Written by Nestor Espinoza, Benjamin Rackham

import glob
from matplotlib.pyplot import figure, plot, subplot, imshow, colorbar, show, title, scatter, text, close, xlim, ylim

# from matplotlib.colors import LogNorm
import numpy as np
import astropy.io.fits as pyfits
import sys

###################################################
# Input Parameters

data_dir = "../data/raw/IMACS"
# dat = '/data/ACCESS/IMACS/WASP103/ut170404_05/'
# file = 'ift0500'
# image = dataset_dir+file
# dataset_dir = dataset_dir
# coords_file = 'w103s_ut170404_05.coords'
markersize = 100

###################################################
# Classes and Functions


def biassec(x):
    x = x.split(",")
    fnum = ((x[0])[1:]).split(":")
    snum = ((x[1])[: len(x[1]) - 1]).split(":")
    fnum[0] = int(float(fnum[0]))
    fnum[1] = int(float(fnum[1]))
    snum[0] = int(float(snum[0]))
    snum[1] = int(float(snum[1]))
    return fnum, snum


def zero_oscan(d):
    zrows = len(np.where(d[:, 0] == 0)[0])
    zcols = len(np.where(d[0, :] == 0)[0])
    if zrows > zcols:
        mode = "r"
        length = d.shape[1]
    else:
        mode = "c"
        length = d.shape[0]
    cmed = []
    for i in range(length):
        if mode == "r":
            data = d[:, i]
        else:
            data = d[i, :]
        I = np.where(data != 0)[0]
        cmed.append(np.median(data[I]))
    return np.median(np.array(cmed))


def BiasTrim(d, c, h, otype, datasec=None):
    """
    Overscan/Bias correct and Trim an IMACS chip
    """
    # bias has no significant structure, so a single median suffices, I think
    # overscan = [0:49] [2097:2145]
    def get_full_datasec(h):
        binning = h["binning"]
        b1, b2 = np.array(binning.split("x")).astype("int")
        datasec = "[1:{:},1:{:}]".format(2048 / b1, 4096 / b2)
        return datasec

    oxsec, oysec = biassec(h["biassec"])
    if datasec == None:
        #        dxsec,dysec = biassec(h['datasec'])
        datasec = get_full_datasec(h)
        dxsec, dysec = biassec(datasec)
    else:
        dxsec, dysec = biassec(datasec)
    if otype == "ift":
        oscan_data = d[(oysec[0] - 1) : oysec[1], (oxsec[0] - 1) : oxsec[1]]
        overscan = np.median(oscan_data)
        if overscan == 0:
            overscan = zero_oscan(oscan_data)
        newdata = d[(dysec[0] - 1) : dysec[1], (dxsec[0] - 1) : dxsec[1]] - overscan
    else:
        d = d.transpose()
        oscan_data = d[oxsec[0] - 1 : oxsec[1], oysec[0] - 1 : oysec[1]]
        overscan = np.median(oscan_data)
        if overscan == 0:
            overscan = zero_oscan(oscan_data)
        newdata = d[dxsec[0] - 1 : dxsec[1], dysec[0] - 1 : dysec[1]] - overscan
    # overscan = np.median(d[:,2048:2112])
    # newdata = d[:4095,0:2048] - overscan
    if (c == "c5") or (c == "c6") or (c == "c7") or (c == "c8"):
        if otype == "iff":
            newdata = newdata[::-1, :]
        else:
            newdata = newdata[::-1, ::-1]

    return newdata


class CoordinateData:
    """
    A simple class to hold coordinate data.
    """

    def __init__(self, filename, skiplines=0):
        self.fname = filename
        self.obj = np.array([])
        self.chip = np.array([])
        self.x = np.array([])
        self.y = np.array([])

        with open(filename) as f:
            for _ in range(skiplines):
                next(f)
            for line in f:
                splitted = line.split()
                if len(splitted) == 0:
                    break
                self.obj = np.append(self.obj, splitted[0])
                self.chip = np.append(self.chip, splitted[1][-1])
                self.x = np.append(self.x, splitted[2])
                self.y = np.append(self.y, splitted[3])
        self.chip = self.chip.astype(int)
        self.x = self.x.astype(float)
        self.y = self.y.astype(float)


###################################################
# Target
targets = glob.glob(data_dir + "*")
string_to_show = "\nChoose target: \n\n"
for i in range(len(targets)):
    string_to_show = string_to_show + "(" + str(i) + ") " + str(targets[i].split("/")[-1]) + "\n"
print(string_to_show)
user_ds = input(">  ")
try:
    target_dir = targets[int(user_ds)] + "/"
except:
    sys.exit("Not a valid target.  Exiting...")

# Dataset
datasets = glob.glob(target_dir + "*ut*")
string_to_show = "\nChoose dataset: \n\n"
for i in range(len(datasets)):
    string_to_show = string_to_show + "(" + str(i) + ") " + str(datasets[i].split("/")[-1]) + "\n"
print(string_to_show)
user_ds = input(">  ")
try:
    dataset_dir = datasets[int(user_ds)] + "/"
except:
    sys.exit("Not a valid dataset.  Exiting...")
files = glob.glob(dataset_dir + "if*c1.fits")
objects = []
string_to_show = "\nChoose file to plot: \n\n"
for i, file in enumerate(files):
    string_to_show += "({:}) FILE: {:} | OBJECT: {:}\n".format(i, file.split("/")[-1], pyfits.getheader(file)["object"])
print(string_to_show)
response = input(">  ")
file = files[int(response)].replace("c1.fits", "")

# Coordinates
coords_files = glob.glob(dataset_dir + "*coords")
string_to_show = "\nChoose coords files: \n\n"
for i in range(len(coords_files)):
    string_to_show = string_to_show + "(" + str(i) + ") " + str(coords_files[i].split("/")[-1]) + "\n"
print(string_to_show)
response = input(">  ")
try:
    coords_file = coords_files[int(response)]
except:
    sys.exit("Not a valid coords file.  Exiting...")

coords = CoordinateData(coords_file)

# Plotting
otype = file.split("/")[-1][0:3]
chips = ["c1", "c2", "c3", "c4", "c6", "c5", "c8", "c7"]
order = [1, 2, 3, 4, 6, 5, 8, 7]
# ranges = [[10,600],[200,500],[300,700],[0,600]]
ranges = [0, 2000]
figure()
for i in range(8):
    #     print 'Plotting chip {:}'.format(i+1)
    subplot(int("24" + str(i + 1)))
    d, h = pyfits.getdata(file + chips[i] + ".fits", header=True)
    if i == 0:
        print("\n Plotting file: {:}".format(file))
        print("    EXPTYPE: {:} | OBJECT:  {:}".format(h["exptype"], h["object"]))
    d = BiasTrim(d, chips[i], h, otype)
    binning = int(h["binning"][0])
    im = imshow(d, interpolation="nearest")  # .transpose())
    for c in range(len(coords.chip)):
        if order[i] == coords.chip[c]:
            scatter(coords.x[c], coords.y[c], s=markersize, marker="x", color="w")
            plot([coords.x[c], coords.x[c]], [0, 4096 / binning], "w", lw=0.25)
            text(coords.x[c], coords.y[c], coords.obj[c], color="w")
    title("Chip " + (chips[i])[1:])
    im.set_clim(ranges[0], ranges[1])
    ylim([d.shape[0], 0])
    xlim([0, d.shape[1]])
show()
close()
