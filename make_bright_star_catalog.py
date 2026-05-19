# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "astropy",
#   "astroquery",
#   "matplotlib",
#   "numpy",
#   "scipy",
# ]
# ///

"""
Downloads Gaia DR3 sources with G < 8.5 and produces a fixed-width text catalog
in TheSkyX user catalog format (sexagesimal RA/Dec, J2000.0).

Filter logic (order matters):
  1. Download ALL Gaia sources with G < 8.5 (no type or variability filter —
     all sources act as contaminants in the proximity checks)
  2. Remove any source that has a brighter source within 2° (field isolation)
  3. Remove any surviving candidate that has ANY other source from the full
     G < 8.5 catalog within 5 arcmin (focus-box isolation)
  4. Final filter: stellar, non-variable, 3.5 <= G <= 6.5
  5. Propagate positions from Gaia J2016.0 to J2000.0 using proper motion
"""

import ssl
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astroquery.gaia import Gaia
import os

# Work around self-signed / intercepted certificates on some networks
ssl._create_default_https_context = ssl._create_unverified_context
import urllib3
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
import requests
_orig_request = requests.Session.request
def _no_verify(self, method, url, **kwargs):
    kwargs.setdefault("verify", False)
    return _orig_request(self, method, url, **kwargs)
requests.Session.request = _no_verify

Gaia.ROW_LIMIT = -1   # remove the default 2000-row cap

OUTPUT_FILE = "bright_star_catalog.txt"
MAG_DOWNLOAD = 8.5    # download limit (2 mag deeper for proximity checks)
MAG_UPPER    = 6.5    # output upper limit
MAG_LOWER    = 3.5    # output lower limit
FIELD_DEG    = 2.0    # field-isolation radius (degrees)
BOX_ARCMIN   = 5.0    # focus-box isolation radius (arcmin)

J2016 = Time("J2016.0")
J2000 = Time("J2000.0")


def query_gaia(mag_limit: float) -> "astropy.table.Table":
    print(f"Querying Gaia DR3 for all sources with G < {mag_limit} ...")
    query = f"""
        SELECT gs.source_id, gs.ra, gs.dec,
               gs.pmra, gs.pmdec,
               gs.phot_g_mean_mag,
               gs.classprob_dsc_combmod_star,
               gs.phot_variable_flag,
               hip.original_ext_source_id AS hip_id
        FROM gaiadr3.gaia_source AS gs
        LEFT JOIN gaiadr3.hipparcos2_best_neighbour AS hip
               ON gs.source_id = hip.source_id
        WHERE gs.phot_g_mean_mag < {mag_limit}
        ORDER BY gs.phot_g_mean_mag
    """
    job = Gaia.launch_job_async(query, verbose=False)
    t = job.get_results()
    print(f"  Retrieved {len(t):,} sources.")
    return t


def skycoords(table):
    return SkyCoord(ra=np.array(table["ra"]) * u.deg,
                    dec=np.array(table["dec"]) * u.deg)


def remove_dominated(candidates, all_sources, sep_deg: float):
    """Remove candidates that have a brighter source within sep_deg degrees."""
    cand_coords = skycoords(candidates)
    all_coords  = skycoords(all_sources)
    cand_mags   = np.array(candidates["phot_g_mean_mag"])
    all_mags    = np.array(all_sources["phot_g_mean_mag"])

    print(f"Field isolation: searching within {sep_deg}° ...")
    idx_cand, idx_all, _, _ = all_coords.search_around_sky(cand_coords, sep_deg * u.deg)

    dominated = set()
    for ic, ia in zip(idx_cand, idx_all):
        if all_mags[ia] < cand_mags[ic]:
            dominated.add(ic)

    keep = np.array([i not in dominated for i in range(len(candidates))])
    print(f"  Removed {int(np.sum(~keep)):,} sources dominated by a brighter neighbour.")
    return candidates[keep]


def remove_crowded(candidates, all_sources, radius_arcmin: float):
    """Remove candidates that have any other source within radius_arcmin arcmin."""
    cand_coords     = skycoords(candidates)
    all_coords      = skycoords(all_sources)
    cand_source_ids = np.array(candidates["source_id"])
    all_ids         = np.array(all_sources["source_id"])

    print(f"Focus-box isolation: searching within {radius_arcmin} arcmin ...")
    radius = radius_arcmin / 60.0 * u.deg
    idx_cand, idx_all, _, _ = all_coords.search_around_sky(cand_coords, radius)

    crowded = set()
    for ic, ia in zip(idx_cand, idx_all):
        if all_ids[ia] != cand_source_ids[ic]:
            crowded.add(ic)

    keep = np.array([i not in crowded for i in range(len(candidates))])
    print(f"  Removed {int(np.sum(~keep)):,} sources with a neighbour within {radius_arcmin} arcmin.")
    return candidates[keep]


def propagate_to_j2000(table) -> tuple[np.ndarray, np.ndarray]:
    """Return (ra_j2000, dec_j2000) arrays propagated from Gaia J2016.0."""
    ra  = np.array(table["ra"],   dtype=float)
    dec = np.array(table["dec"],  dtype=float)
    pmra  = np.array(table["pmra"],  dtype=float)
    pmdec = np.array(table["pmdec"], dtype=float)

    # Stars with no proper motion measurement — treat as zero
    pmra  = np.where(np.isfinite(pmra),  pmra,  0.0)
    pmdec = np.where(np.isfinite(pmdec), pmdec, 0.0)

    coords = SkyCoord(
        ra=ra * u.deg,
        dec=dec * u.deg,
        pm_ra_cosdec=pmra  * u.mas / u.yr,
        pm_dec=pmdec * u.mas / u.yr,
        obstime=J2016,
        frame="icrs",
    )
    import warnings
    from erfa import ErfaWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        coords_j2000 = coords.apply_space_motion(new_obstime=J2000)
    return coords_j2000.ra.deg, coords_j2000.dec.deg


def to_hms(ra_deg):
    total_s = ra_deg * 3600.0 / 15.0
    h = int(total_s // 3600)
    total_s -= h * 3600
    m = int(total_s // 60)
    s = total_s - m * 60
    return h, m, s


def to_dms(dec_deg):
    sign = "+" if dec_deg >= 0 else "-"
    total_s = abs(dec_deg) * 3600.0
    d = int(total_s // 3600)
    total_s -= d * 3600
    m = int(total_s // 60)
    s = total_s - m * 60
    return sign, d, m, s


def format_name(row) -> str:
    hip = row["hip_id"]
    if hip is not None and str(hip).strip() not in ("", "--", "nan"):
        return f"HIP {int(hip)}"
    return f"GAI {row['source_id']}"


def plot_sky(table, ra_j2000, dec_j2000, path: str = "bright_star_catalog.png"):
    import matplotlib.pyplot as plt

    mags = np.array(table["phot_g_mean_mag"], dtype=float)

    # Mollweide expects RA in radians, wrapped to [-π, π] (east left)
    ra_rad  = np.deg2rad(ra_j2000 - 180.0)
    dec_rad = np.deg2rad(dec_j2000)

    sizes = 60 * (MAG_UPPER - mags + 1) ** 1.5

    fig = plt.figure(figsize=(14, 7))
    ax  = fig.add_subplot(111, projection="mollweide")
    ax.scatter(ra_rad, dec_rad, s=sizes, c=mags, cmap="plasma_r",
               alpha=0.7, linewidths=0, vmin=MAG_LOWER, vmax=MAG_UPPER)

    cb = plt.colorbar(ax.collections[0], ax=ax, orientation="horizontal",
                      pad=0.05, fraction=0.03)
    cb.set_label("Gaia G magnitude")

    ax.set_title(f"Focus star catalog  ({len(table):,} stars,  "
                 f"{MAG_LOWER} ≤ G ≤ {MAG_UPPER},  J2000.0)", pad=16)
    ax.grid(True, alpha=0.3)
    ax.set_xticklabels(["22h", "20h", "18h", "16h", "14h",
                        "12h", "10h",  "8h",  "6h",  "4h",  "2h"])

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    print(f"Sky plot saved to: {os.path.abspath(path)}")
    plt.show()


def write_catalog(table, ra_j2000, dec_j2000, path: str):
    W_RAH, W_RAM, W_RAS = 4, 4, 7
    W_SGN, W_DED, W_DEM, W_DES = 2, 4, 4, 6
    W_MAG  = 7
    W_TYPE = 10
    W_NAME = 16

    header = " ".join([
        f"{'RAh':>{W_RAH}}", f"{'RAm':>{W_RAM}}", f"{'RAs':>{W_RAS}}",
        f"{'Sgn':>{W_SGN}}", f"{'DEd':>{W_DED}}", f"{'DEm':>{W_DEM}}",
        f"{'DEs':>{W_DES}}", f"{'Mag':>{W_MAG}}", f"{'ObjType':>{W_TYPE}}",
        f"{'Label/Search':>{W_NAME}}",
    ])
    separator = "-" * len(header)

    print(f"Writing {len(table):,} stars to {path} ...")
    with open(path, "w") as fh:
        fh.write(header + "\n")
        fh.write(separator + "\n")
        for i, row in enumerate(table):
            h, m, s       = to_hms(ra_j2000[i])
            sg, d, dm, ds = to_dms(dec_j2000[i])
            mag  = float(row["phot_g_mean_mag"])
            name = format_name(row)
            line = " ".join([
                f"{h:>{W_RAH}d}", f"{m:>{W_RAM}d}", f"{s:>{W_RAS}.3f}",
                f"{sg:>{W_SGN}}", f"{d:>{W_DED}d}", f"{dm:>{W_DEM}d}",
                f"{ds:>{W_DES}.2f}", f"{mag:>{W_MAG}.3f}",
                f"{'Star':>{W_TYPE}}", f"{name:>{W_NAME}}",
            ])
            fh.write(line + "\n")

    print(f"Done. Catalog written to: {os.path.abspath(path)}")


def main():
    # Step 1 — download all sources to MAG_DOWNLOAD
    all_sources = query_gaia(MAG_DOWNLOAD)

    # Step 2 — field isolation
    candidates = remove_dominated(all_sources, all_sources, FIELD_DEG)

    # Step 3 — focus-box isolation (checked against full download)
    candidates = remove_crowded(candidates, all_sources, BOX_ARCMIN)

    # Step 4 — final filter: stellar, non-variable, within magnitude range
    mags     = np.array(candidates["phot_g_mean_mag"])
    probs    = np.array(candidates["classprob_dsc_combmod_star"])
    var_flag = np.array(candidates["phot_variable_flag"], dtype=str)

    is_star    = (probs > 0.5) | np.isnan(probs)
    not_var    = var_flag != "VARIABLE"
    in_range   = (mags >= MAG_LOWER) & (mags <= MAG_UPPER)
    candidates = candidates[is_star & not_var & in_range]
    print(f"After final filter (stellar, non-variable, "
          f"{MAG_LOWER} <= G <= {MAG_UPPER}): {len(candidates):,} stars remain.")

    # Step 5 — propagate positions to J2000.0
    print("Propagating positions from J2016.0 to J2000.0 ...")
    ra_j2000, dec_j2000 = propagate_to_j2000(candidates)

    write_catalog(candidates, ra_j2000, dec_j2000, OUTPUT_FILE)
    plot_sky(candidates, ra_j2000, dec_j2000)


if __name__ == "__main__":
    main()
