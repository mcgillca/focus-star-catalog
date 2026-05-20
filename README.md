# Focus Star Catalog

A catalog of ~793 isolated bright stars suitable for telescope focusing, generated from Gaia DR3. Designed for use with TheSkyX as a user-added catalog.

Optimised for the Takahashi FSQ-85ED with a Player One Zeus (IMX455) camera, which gives a 4.58° × 3.05° field of view (5.5° diagonal). The 3° field isolation radius guarantees no brighter star appears anywhere in the field, including the corners.

![Sky distribution of focus stars](bright_star_catalog.png)

## Catalog properties

- **Stars**: ~793 all-sky
- **Magnitude range**: 3.5 ≤ G ≤ 6.5
- **Field isolation**: no brighter Gaia source within 3° — covers the full FOV diagonal of the FSQ-85/Zeus combination
- **Focus-box isolation**: no other Gaia source (to G < 9.5, ~16× fainter than the brightest focus star) within 5 arcmin — clean PSF for HFD/centroid algorithms
- **Non-variable**: Gaia `phot_variable_flag != VARIABLE`
- **Stellar sources only**: `classprob_dsc_combmod_star > 0.5` (or NULL for bright stars not assessed by the pipeline)
- **Coordinates**: sexagesimal (H M S / sign D M S), ICRS J2000.0, propagated from Gaia DR3 J2016.0 using per-star proper motions
- **Names**: HIP numbers where available (Hipparcos cross-match), Gaia source_id otherwise

## Files

| File | Description |
|------|-------------|
| `make_bright_star_catalog.py` | Script to regenerate the catalog from Gaia DR3 |
| `bright_star_catalog.txt` | Catalog in TheSkyX fixed-width text format (sexagesimal RA/Dec, J2000.0) |
| `bright_star_catalog.SDBX` | Pre-built TheSkyX user catalog database |
| `bright_star_catalog.png` | Mollweide all-sky plot of catalog stars |
| `AtFocus2.dbq` | Updated database query file used by @Focus2 to select focus stars |

## Regenerating the catalog

Requires [uv](https://docs.astral.sh/uv/). Dependencies are managed automatically via PEP 723 inline metadata.

```bash
uv run make_bright_star_catalog.py
```

This will query the Gaia DR3 archive (~30–60 seconds for the async job), apply the isolation filters, and write a fresh `bright_star_catalog.txt` and `bright_star_catalog.png`.

## Importing into TheSkyX

1. **Back up the original database query file** before making any changes:
   ```
   ~/Library/Application Support/Software Bisque/TheSkyX Professional Edition/Database Queries/AtFocus2.dbq
   ```

2. Copy `AtFocus2.dbq` to:
   ```
   ~/Library/Application Support/Software Bisque/TheSkyX Professional Edition/Database Queries/
   ```

3. Copy `bright_star_catalog.SDBX` to:
   ```
   ~/Library/Application Support/Software Bisque/TheSkyX Professional Edition/SDBs/
   ```

If you have regenerated the catalog from the script, import `bright_star_catalog.txt` into TheSkyX to produce a new `.SDBX` file, then copy that into the `SDBs` folder.

Stars will then appear in TheSkyX with an `AF2` prefix. When searching, include the prefix — for example:

| Star | Search for |
|------|-----------|
| HIP 36046 | `AF2 HIP 36046` |
| Gaia source 3314024566919613952 | `AF2 GAI 3314024566919613952` |

## Customising the catalog

The key parameters are on lines 48–52 of `make_bright_star_catalog.py`:

```python
MAG_DOWNLOAD = 9.5    # download limit (3 mag deeper, ~16x fainter than brightest focus star)
MAG_UPPER    = 6.5    # output upper limit
MAG_LOWER    = 3.5    # output lower limit
FIELD_DEG    = 3.0    # field-isolation radius (degrees) — covers FSQ-85/Zeus FOV corners
BOX_ARCMIN   = 5.0    # focus-box isolation radius (arcmin)
```

After adjusting these and running the script, a new `bright_star_catalog.txt` will be produced. Import it into TheSkyX to generate a new `.SDBX` database file, then copy that into the `SDBs` folder as described above.

## Filter logic

The script applies filters in a specific order to ensure sources outside the final magnitude range still act as contaminants during proximity checks:

1. Download all Gaia sources with G < 9.5 (no type or variability filtering at this stage)
2. Propagate all positions from Gaia J2016.0 to today using per-star proper motions
3. Remove any source with a brighter source within 3° (field isolation, using today's positions)
4. Remove any surviving candidate with any other source within 5 arcmin (focus-box isolation, checked against the full G < 9.5 download)
5. Keep only stellar, non-variable sources with 3.5 ≤ G ≤ 6.5
6. Fetch Hipparcos cross-matches for the final stars
7. Propagate output positions from Gaia J2016.0 to J2000.0 for the catalog
