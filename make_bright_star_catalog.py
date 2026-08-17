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
Downloads Gaia DR3 sources with G < 9.5 and produces a fixed-width text catalog
in TheSkyX user catalog format (sexagesimal RA/Dec, CATALOG_EPOCH).

Filter logic (order matters):
  1. Download ALL Gaia sources with G < 9.5 (no type or variability filter —
     all sources act as contaminants in the proximity checks)
  2. Propagate all positions from Gaia J2016.0 to CATALOG_EPOCH using proper motions
  3. Remove any source that has a brighter source within 3° (field isolation)
  4. Remove any surviving candidate that has ANY other source from the full
     G < 9.5 catalog within 5 arcmin (focus-box isolation)
  5. Final filter: stellar, non-variable, 3.5 <= G <= 6.5
  6. Output positions are at CATALOG_EPOCH; TheSkyX handles precession internally
"""

import ssl
import warnings
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

OUTPUT_FILE  = "bright_star_catalog.txt"
MAG_DOWNLOAD = 9.5    # download limit (3 mag deeper, ~16x fainter than brightest focus star)
MAG_UPPER    = 6.5    # output upper limit
MAG_LOWER    = 3.5    # output lower limit
FIELD_DEG    = 3.0    # field-isolation radius (degrees) — covers FSQ-85/Zeus FOV corners
BOX_ARCMIN   = 5.0    # focus-box isolation radius (arcmin)

T_GAIA         = Time("J2016.0")   # Gaia DR3 reference epoch
CATALOG_EPOCH  = Time("J2030.0")   # target epoch for filtering and output


def query_gaia(mag_limit: float) -> "astropy.table.Table":
    """Download all sources for proximity filtering (no HIP join)."""
    print(f"Querying Gaia DR3 for all sources with G < {mag_limit} ...")
    query = f"""
        SELECT source_id, ra, dec, pmra, pmdec,
               phot_g_mean_mag,
               classprob_dsc_combmod_star,
               phot_variable_flag
        FROM gaiadr3.gaia_source
        WHERE phot_g_mean_mag < {mag_limit}
        ORDER BY phot_g_mean_mag
    """
    job = Gaia.launch_job_async(query, verbose=False)
    t = job.get_results()
    print(f"  Retrieved {len(t):,} sources.")
    return t


def query_hip_ids(source_ids: np.ndarray) -> dict:
    """Fetch Hipparcos IDs for a small list of source_ids."""
    print(f"Fetching HIP cross-matches for {len(source_ids):,} final stars ...")
    id_list = ", ".join(str(s) for s in source_ids)
    query = f"""
        SELECT source_id, original_ext_source_id AS hip_id
        FROM gaiadr3.hipparcos2_best_neighbour
        WHERE source_id IN ({id_list})
    """
    job = Gaia.launch_job(query, verbose=False)
    t = job.get_results()
    return {int(row["source_id"]): int(row["hip_id"]) for row in t}


def propagate(ra_deg, dec_deg, pmra_masyr, pmdec_masyr,
              t_from: Time, t_to: Time) -> tuple[np.ndarray, np.ndarray]:
    """Propagate positions from t_from to t_to using proper motions."""
    pmra_masyr  = np.where(np.isfinite(pmra_masyr),  pmra_masyr,  0.0)
    pmdec_masyr = np.where(np.isfinite(pmdec_masyr), pmdec_masyr, 0.0)

    coords = SkyCoord(
        ra=ra_deg * u.deg,
        dec=dec_deg * u.deg,
        pm_ra_cosdec=pmra_masyr  * u.mas / u.yr,
        pm_dec=pmdec_masyr * u.mas / u.yr,
        obstime=t_from,
        frame="icrs",
    )
    from erfa import ErfaWarning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ErfaWarning)
        coords_new = coords.apply_space_motion(new_obstime=t_to)
    return coords_new.ra.deg, coords_new.dec.deg


def skycoords_from_arrays(ra_deg, dec_deg) -> SkyCoord:
    return SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg)


def remove_dominated(cand_ra, cand_dec, cand_mags,
                     all_ra, all_dec, all_mags, sep_deg: float) -> np.ndarray:
    """Return boolean keep-mask: remove candidates with a brighter source within sep_deg."""
    cand_coords = skycoords_from_arrays(cand_ra, cand_dec)
    all_coords  = skycoords_from_arrays(all_ra,  all_dec)

    print(f"Field isolation: searching within {sep_deg}° ...")
    idx_cand, idx_all, _, _ = all_coords.search_around_sky(cand_coords, sep_deg * u.deg)

    dominated = set()
    for ic, ia in zip(idx_cand, idx_all):
        if all_mags[ia] < cand_mags[ic]:
            dominated.add(ic)

    keep = np.array([i not in dominated for i in range(len(cand_ra))])
    print(f"  Removed {int(np.sum(~keep)):,} sources dominated by a brighter neighbour.")
    return keep


def remove_crowded(cand_ra, cand_dec, cand_ids,
                   all_ra, all_dec, all_ids, radius_arcmin: float) -> np.ndarray:
    """Return boolean keep-mask: remove candidates with any other source within radius_arcmin."""
    cand_coords = skycoords_from_arrays(cand_ra, cand_dec)
    all_coords  = skycoords_from_arrays(all_ra,  all_dec)

    print(f"Focus-box isolation: searching within {radius_arcmin} arcmin ...")
    radius = radius_arcmin / 60.0 * u.deg
    idx_cand, idx_all, _, _ = all_coords.search_around_sky(cand_coords, radius)

    crowded = set()
    for ic, ia in zip(idx_cand, idx_all):
        if all_ids[ia] != cand_ids[ic]:
            crowded.add(ic)

    keep = np.array([i not in crowded for i in range(len(cand_ra))])
    print(f"  Removed {int(np.sum(~keep)):,} sources with a neighbour within {radius_arcmin} arcmin.")
    return keep


def format_name(source_id: int, hip_map: dict) -> str:
    hip = hip_map.get(source_id)
    if hip is not None:
        return f"HIP {hip}"
    return f"GAI {source_id}"


def plot_sky(table, ra, dec, path: str = "bright_star_catalog.png"):
    import matplotlib.pyplot as plt

    mags = np.array(table["phot_g_mean_mag"], dtype=float)

    # Mollweide expects RA in radians, wrapped to [-π, π] (east left)
    ra_rad  = np.deg2rad(ra - 180.0)
    dec_rad = np.deg2rad(dec)
    sizes   = 60 * (MAG_UPPER - mags + 1) ** 1.5

    fig = plt.figure(figsize=(14, 7))
    ax  = fig.add_subplot(111, projection="mollweide")
    ax.scatter(ra_rad, dec_rad, s=sizes, c=mags, cmap="plasma_r",
               alpha=0.7, linewidths=0, vmin=MAG_LOWER, vmax=MAG_UPPER)

    cb = plt.colorbar(ax.collections[0], ax=ax, orientation="horizontal",
                      pad=0.05, fraction=0.03)
    cb.set_label("Gaia G magnitude")

    ax.set_title(f"Focus star catalog  ({len(table):,} stars,  "
                 f"{MAG_LOWER} ≤ G ≤ {MAG_UPPER},  {CATALOG_EPOCH.value})", pad=16)
    ax.grid(True, alpha=0.3)
    ax.set_xticklabels(["22h", "20h", "18h", "16h", "14h",
                        "12h", "10h",  "8h",  "6h",  "4h",  "2h"])

    plt.tight_layout()
    plt.savefig(path, dpi=150)
    print(f"Sky plot saved to: {os.path.abspath(path)}")
    plt.show()


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


def write_catalog(table, ra, dec, hip_map: dict, path: str):
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
            h, m, s       = to_hms(ra[i])
            sg, d, dm, ds = to_dms(dec[i])
            mag  = float(row["phot_g_mean_mag"])
            name = format_name(int(row["source_id"]), hip_map)
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

    # Step 2 — propagate ALL positions from J2016.0 to CATALOG_EPOCH
    print(f"Propagating positions from J2016.0 to {CATALOG_EPOCH.value} ...")
    all_ra, all_dec = propagate(
        np.array(all_sources["ra"],    dtype=float),
        np.array(all_sources["dec"],   dtype=float),
        np.array(all_sources["pmra"],  dtype=float),
        np.array(all_sources["pmdec"], dtype=float),
        T_GAIA, CATALOG_EPOCH,
    )
    all_mags = np.array(all_sources["phot_g_mean_mag"], dtype=float)
    all_ids  = np.array(all_sources["source_id"])

    # Step 3 — field isolation using propagated positions
    keep = remove_dominated(all_ra, all_dec, all_mags,
                            all_ra, all_dec, all_mags, FIELD_DEG)
    cand_ra, cand_dec = all_ra[keep], all_dec[keep]
    cand_mags = all_mags[keep]
    cand_ids  = all_ids[keep]
    candidates = all_sources[keep]

    # Step 4 — focus-box isolation (checked against full propagated download)
    keep = remove_crowded(cand_ra, cand_dec, cand_ids,
                          all_ra,  all_dec,  all_ids, BOX_ARCMIN)
    cand_ra, cand_dec = cand_ra[keep], cand_dec[keep]
    candidates = candidates[keep]

    # Step 5 — final filter: stellar, non-variable, within magnitude range
    mags     = np.array(candidates["phot_g_mean_mag"], dtype=float)
    probs    = np.array(candidates["classprob_dsc_combmod_star"], dtype=float)
    var_flag = np.array(candidates["phot_variable_flag"], dtype=str)

    is_star  = (probs > 0.5) | np.isnan(probs)
    not_var  = var_flag != "VARIABLE"
    in_range = (mags >= MAG_LOWER) & (mags <= MAG_UPPER)
    mask     = is_star & not_var & in_range

    candidates = candidates[mask]
    cand_ra    = cand_ra[mask]
    cand_dec   = cand_dec[mask]
    print(f"After final filter (stellar, non-variable, "
          f"{MAG_LOWER} <= G <= {MAG_UPPER}): {len(candidates):,} stars remain.")

    # Step 6 — fetch HIP IDs for the final stars only (small targeted query)
    hip_map = query_hip_ids(np.array(candidates["source_id"]))

    # Step 7 — output uses the already-propagated CATALOG_EPOCH positions
    write_catalog(candidates, cand_ra, cand_dec, hip_map, OUTPUT_FILE)
    plot_sky(candidates, cand_ra, cand_dec)


if __name__ == "__main__":
    main()
