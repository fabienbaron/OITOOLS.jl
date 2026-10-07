# PMOIRED example models

Small models in PMOIRED's dict syntax, for trying the Model perspective's **Import PMOIRED…**
without having to write one first.

| file | model |
|---|---|
| `uniform_disc.py` | a single uniform disc |
| `limb_darkened_power_law.py` | a power-law limb-darkened disc, `I(μ) = μ^α` |
| `binary.py` | two uniform discs with a separation and a flux ratio |
| `star_plus_gaussian_ring.py` | an unresolved star inside an inclined Gaussian ring |

Each was written by `dict_to_pmoired_file` and read back with `pmoired_to_dict`, so every one
is known to import. That matters more than it sounds: the importer is a **transpiler, not a
validator**, so a file it cannot read fails at the point of use rather than being rejected —
and an example that does not import is worse than no example at all.

Two conventions differ between the packages and are worth knowing before importing anything of
your own:

* azimuthal modes — OITOOLS uses `+π/2` where PMOIRED uses `−π/2`, so `az projang` needs
  adjusting by hand;
* the OITOOLS-only geometries (`ldlin`, `ldquad`, `ldpow`, `resolved`) have no PMOIRED
  equivalent, and exporting one warns rather than silently writing something that means
  something else on the other side.

After importing, look at the Creation page's inspector: an unrecognised key does not error, it
silently changes what the component IS.
