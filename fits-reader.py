from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import ndimage
import numpy as np
import math

# Reading FITS files
hdu = fits.open('2014/kanz_bband_fi_20140226_114908.fts')
data = hdu[0].data
header = hdu[0].header


# Visualization FITS data

'''
For inversing color in colorbar, add this:

orig_map = plt.cm.get_cmap('gray')
reversed_map = orig_map.reversed()
'''
# plt.figure(figsize=(13, 10))
# plt.imshow(data, cmap='gray', origin='lower')
# plt.colorbar()
# plt.show()

# Erosion and filtering data

X = [i for i in range(data.shape[1])]
filter_data = ndimage.grey_erosion(data, size=(13, 13))
sub_data = data - filter_data
I = np.max(sub_data)/2
sub_data = np.where(sub_data > I, sub_data, 0)
for i in range(data.shape[0]):
    for j in range(data.shape[1]):
        if j < 400 or j > 1600:
            sub_data[i][j] = 0
for i in range(data.shape[0]):
    if i < 400 or i > 1600:
        for j in range(data.shape[1]):
                sub_data[i][j] = 0

y, x = np.nonzero(sub_data)

# Defining the latitude of the Sunspots
x_0 = header['CRPIX1'] # X coord of center of the Sun disk in px
y_0 = header['CRPIX2'] # Y coord of center of the Sun disk in px
dx = 4.8e-6 * header['CDELT1'] # Angular size along X in rad
dy = 4.8e-6 * header['CDELT2'] # Angular size along Y in rad
P0 = header['SOLAR_P0']/180*math.pi # Angle of the Sun rotational axis
B0 = header['SOLAR_B0']/180 * math.pi # Latitude of the center of the Sun disk in rad
L0 = header['ANGLE']/180 * math.pi # Longitude of the center of the Sun disk in rad
alpha_0 = header['CRVAL1'] * 4.8e-6 # RA of center of the Sun disk in rad
delta_0 = header['CRVAL2'] * 4.8e-6 # Declination of center of the Sun disk in rad

def thetax_thetay(x, y):
    X = dx * x * math.cos(P0) - dy * y * math.sin(P0)
    Y = dx * x * math.sin(P0) + dy * y * math.cos(P0)
    thetax = np.arctan(math.pi * X / 180)
    thetay = np.arctan(math.pi/180 * Y / np.sqrt(1 + (math.pi*X/180)**2))
    return thetax, thetay
def rho_phi(x, y):
    X = dx * x * math.cos(P0) - dy * y * math.sin(P0)
    Y = dx * x * math.sin(P0) + dy * y * math.cos(P0)
    tanx =  X
    tany =  Y / np.sqrt(1 + X**2)
    phi = P0 + np.arctan(-X/Y)  # Position angle
    rho = np.sqrt(X**2 + Y**2)
    return rho, phi

def carrington_heliogr_coord(x, y):
    rho, phi = rho_phi(x, y)
    B = np.arcsin(math.sin(B0) * np.cos(rho) + math.cos(B0) * np.sin(rho) * np.cos(P0 - phi))
    L = np.arcsin(np.sin(rho) * np.sin(P0 - phi) / math.cos(B0)) + L0
    return B * 180 / math.pi, L * 180 / math.pi

B, L = carrington_heliogr_coord(x, y)

print('y[i]_px', '\t', 'x[i]_px', '\t', 'B[i]_deg', '\t', 'L[i]_deg')
for i in range(len(B)):
    print(y[i], '\t', x[i], '\t', B[i], '\t', L[i])

# Plotting Sunspots from x and y
f1 = plt.figure(figsize=[13,10])
im1 = plt.imshow(data, cmap='gray', origin='lower')
plt.colorbar(im1)

im0 = plt.scatter(x, y, s = 10, c=[sub_data[y[i]][x[i]] for i in range(len(x))],cmap='jet')
plt.colorbar(im0)
plt.xlabel(r'X-axis, px')
plt.ylabel(r'Y-axis, px')
plt.title(r'Sunspots with P0=' + f'{P0.__round__(2)}' + '; B0=' + f'{B0.__round__(2)}' + '; L0=' + f'{L0.__round__(2)}',fontsize=11)
plt.show()

# Plotting all operations with the FITS data

for k in range(659, 1600):
    fig = plt.figure(figsize=(10, 7))
    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)
    # Original datass
    fig.suptitle(f'Number of raw is {k}')
    ax1.plot(X, -data[k][:])
    # Erosion data
    ax2.plot(X, -filter_data[k][:])
    # Combined image
    ax3.plot(X, -data[k][:], -filter_data[k][:])
    # Subtracted data
    ax4.plot(X, sub_data[k][:])
    plt.show(block=True)
    plt.pause(0.5)
    plt.close()
