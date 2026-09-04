#!/usr/bin/env python3
"""Compute chi2 of 2004true64.fits against 2004-data1.oifits using PMOIRED's visImage logic."""
import numpy as np
from astropy.io import fits as pyfits
from astropy.io.fits import open as fitsopen

# === PMOIRED's visImage (inlined from oifake.py:801) ===
def visImage(image, scale, u, v, wl):
    """
    image: 2D image, square pixels
    scale: pixel size in mas
    u, v: spatial coordinates (vectors, in m)
    wl: wavelength (scalar or vector, in um)
    """
    Ny, Nx = np.shape(image)
    x = np.arange(Nx) - (Nx-1)/2
    y = np.arange(Ny) - (Ny-1)/2
    X, Y = np.meshgrid(scale*x, scale*y)

    c = np.pi / 180 / 3600 / 1000 / 1e-6
    if np.isscalar(wl):
        vis = image[:, :, None] * \
            np.exp(-2j * np.pi * c * (X[:, :, None] * u[None, None, :] +
                                       Y[:, :, None] * v[None, None, :]) / wl)
    else:
        vis = image[:, :, None, None] * \
            np.exp(-2j * np.pi * c * (X[:, :, None, None] * u[None, None, :, None] +
                                       Y[:, :, None, None] * v[None, None, :, None]) /
                   wl[None, None, None, :])
    return np.sum(vis, axis=(0, 1)) / np.sum(image)

# === Load OIFITS using astropy ===
def load_oifits(filename):
    """Minimal OIFITS reader - extract V2 and T3 data."""
    hdul = pyfits.open(filename)
    wl_um = None
    v2_list = []
    t3_list = []

    for hdu in hdul:
        if hdu.name == 'OI_WAVELENGTH':
            wl_um = hdu.data['EFF_WAVE'] * 1e6  # m -> um
        elif hdu.name == 'OI_VIS2':
            d = hdu.data
            for row in d:
                v2_list.append({
                    'V2': row['VIS2DATA'], 'EV2': row['VIS2ERR'],
                    'u': row['UCOORD'], 'v': row['VCOORD'],
                    'flag': row['FLAG'] if 'FLAG' in d.columns.names else False
                })
        elif hdu.name == 'OI_T3':
            d = hdu.data
            for row in d:
                t3_list.append({
                    'T3PHI': row['T3PHI'], 'T3PHIERR': row['T3PHIERR'],
                    'T3AMP': row['T3AMP'], 'T3AMPERR': row['T3AMPERR'],
                    'u1': row['U1COORD'], 'v1': row['V1COORD'],
                    'u2': row['U2COORD'], 'v2': row['V2COORD'],
                    'flag': row['FLAG'] if 'FLAG' in d.columns.names else False
                })
    hdul.close()
    return wl_um, v2_list, t3_list

# === Main ===
# Load image
hdul = pyfits.open('data/2004true64.fits')
image = hdul[0].data.astype(np.float64)
header = hdul[0].header
hdul.close()
print(f"Image shape: {image.shape}, sum={image.sum():.6f}")
print(f"CDELT1={header.get('CDELT1')}, CUNIT1={header.get('CUNIT1')}")

scale = abs(header['CDELT1']) * 1000  # arcsec -> mas
print(f"Pixel scale: {scale} mas")

image_norm = image / image.sum()

# Load OIFITS
wl_um, v2_list, t3_list = load_oifits('data/BC2004/2004-data1.oifits')
print(f"Wavelength: {wl_um} um")
print(f"N_V2={len(v2_list)}, N_T3={len(t3_list)}")

# For each wavelength channel
wl_scalar = wl_um[0]  # monochromatic

# === V2 chi2 ===
chi2_v2 = 0.0
nv2 = 0
for d in v2_list:
    flag = d['flag']
    if np.any(flag):
        continue
    u_arr = np.array([d['u']])
    v_arr = np.array([d['v']])
    cvis = visImage(image_norm, scale, u_arr, v_arr, wl_scalar)
    V2_model = float(np.abs(cvis[0])**2)
    V2_data = float(np.atleast_1d(d['V2'])[0])
    EV2 = float(np.atleast_1d(d['EV2'])[0])

    if EV2 > 0:
        resid = (V2_model - V2_data) / EV2
        chi2_v2 += resid**2
        nv2 += 1
        if nv2 <= 5:
            print(f"  V2[{nv2}]: data={V2_data:.6f} model={V2_model:.6f} err={EV2:.6f} resid={resid:.2f}")

print(f"\nV2: chi2={chi2_v2:.2f}  nv2={nv2}  chi2/nv2={chi2_v2/nv2:.2f}")

# === T3phi chi2 ===
chi2_t3phi = 0.0
nt3 = 0
for d in t3_list:
    flag = d['flag']
    if np.any(flag):
        continue
    u1, v1 = np.array([d['u1']]), np.array([d['v1']])
    u2, v2 = np.array([d['u2']]), np.array([d['v2']])
    u3, v3 = -(u1 + u2), -(v1 + v2)

    cvis1 = visImage(image_norm, scale, u1, v1, wl_scalar)
    cvis2 = visImage(image_norm, scale, u2, v2, wl_scalar)
    cvis3 = visImage(image_norm, scale, u3, v3, wl_scalar)

    t3 = cvis1[0] * cvis2[0] * cvis3[0]
    T3phi_model = float(np.angle(t3) * 180 / np.pi)
    T3phi_data = float(np.atleast_1d(d['T3PHI'])[0])
    ET3phi = float(np.atleast_1d(d['T3PHIERR'])[0])

    if ET3phi > 0:
        diff = T3phi_model - T3phi_data
        diff = (diff + 180) % 360 - 180
        resid = diff / ET3phi
        chi2_t3phi += resid**2
        nt3 += 1
        if nt3 <= 5:
            print(f"  T3phi[{nt3}]: data={T3phi_data:.2f} model={T3phi_model:.2f} err={ET3phi:.2f} resid={resid:.2f}")

print(f"\nT3phi: chi2={chi2_t3phi:.2f}  nt3={nt3}  chi2/nt3={chi2_t3phi/nt3:.2f}")
print(f"\nTotal: chi2={chi2_v2+chi2_t3phi:.2f}  N={nv2+nt3}  chi2/N={(chi2_v2+chi2_t3phi)/(nv2+nt3):.2f}")
