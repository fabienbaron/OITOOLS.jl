#!/usr/bin/env bash
#
# Package the built bundle as a single-file AppImage.
#
#     julia --project=bin app/build.jl      # first: produce build/OITOOLS
#     app/appimage.sh                       # then: OITOOLS-x86_64.AppImage
#
# An AppImage is a squashfs image with a small runtime prepended: running it mounts the image
# and executes AppRun inside. So the user gets ONE executable file rather than the 2.6 GB
# directory `create_app` produces, and nothing is installed or unpacked.
#
# WHY THE PATHS STILL WORK. The mount point is a fresh /tmp/.mount_XXXXXX every launch, so
# nothing may assume a fixed location -- which is exactly the property `resource()` was built
# for: it anchors on `<Sys.BINDIR>/../share/oitools`, and Sys.BINDIR follows the mount. The
# same is true of the Makie font atlas, which is seeded into the user's own cache rather than
# written beside the binary.
#
# COMPRESSION IS A TRADE, not a free win. The image is decompressed BLOCK BY BLOCK as it is
# read, so a smaller file costs CPU on every launch -- and this payload is a 1.2 GB sysimage
# that gets mapped at startup.
#
# THE COMPRESSOR IS NOT A CHOICE HERE. The appimagetool this script downloads carries its own
# mksquashfs, and that build supports **zstd only** -- `--comp xz` fails outright with
# "Compressor \"xz\" is not supported". OITOOLS_APPIMAGE_COMP exists for a host mksquashfs that
# has more, and is otherwise nothing to turn.
#
# THE TWO LEVERS THAT DO WORK ARE THE BLOCK SIZE AND THE LEVEL, AND THEY ARE NOT INDEPENDENT.
# squashfs compresses each block on its own, so at the 128K default a high zstd level has too
# little data in front of it to find matches in. Measured on this bundle (2.6 GB staged):
#
#     blocks   level     AppImage    packing
#     128K     15        589.0 MB     11 s     <- appimagetool's defaults
#     1M       15        569.3 MB      9 s     -19.7 MB
#     128K     19        582.7 MB     27 s      -6.3 MB
#     1M       19        523.2 MB     25 s     -65.8 MB, and the two below add to -26
#
# So neither alone is worth much and together they take 11.2% off: the block size is what lets
# the level pay. Both are defaulted on below. Level costs build time and NOT launch time --
# zstd decompression speed barely depends on it -- and a 1M block is read-ahead rather than
# waste here, because a sysimage is read broadly rather than in scattered fragments.
#
# Not a compression trick, but the real bandwidth win for repeat users: `-u <update string>`
# writes a zsync file, and a later AppImage is then fetched as a DELTA against the one already
# on disk. It needs a published URL, so it belongs to whatever hosts the releases.

set -euo pipefail

here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
root=$(dirname "$here")
bundle=${OITOOLS_APP_DIR:-$root/build/OITOOLS}
comp=${OITOOLS_APPIMAGE_COMP:-zstd}
block=${OITOOLS_APPIMAGE_BLOCK:-1M}
level=${OITOOLS_APPIMAGE_LEVEL:-19}
out=${OITOOLS_APPIMAGE_OUT:-$root/build/OITOOLS-x86_64.AppImage}
cache=${XDG_CACHE_HOME:-$HOME/.cache}/oitools-build

[ -x "$bundle/OITOOLS" ] || { echo "no bundle at $bundle -- run app/build.jl first" >&2; exit 1; }
mkdir -p "$cache"

# ── the tool ─────────────────────────────────────────────────────────────────
tool=$cache/appimagetool-x86_64.AppImage
if [ ! -x "$tool" ]; then
    echo "fetching appimagetool..."
    curl -fsSL -o "$tool" \
        https://github.com/AppImage/appimagetool/releases/download/continuous/appimagetool-x86_64.AppImage
    chmod +x "$tool"
fi

# ── the AppDir ───────────────────────────────────────────────────────────────
#
# AppRun, one .desktop and one icon must sit at the top; everything else is free-form, so the
# bundle goes in whole. The icon's basename must match the desktop file's Icon= key.
appdir=$cache/OITOOLS.AppDir
rm -rf "$appdir"; mkdir -p "$appdir"
cp -a "$bundle"/. "$appdir"/
cp "$root/app/assets/oitools.png" "$appdir/oitools.png"

# AppRun defers to the bundle's own launcher, which is what decides the windowing system --
# GLFW and Qt are both initialised from the sysimage before julia_main runs, so it cannot be
# chosen later. Invoking it by full path makes its own `dirname $0` resolve to $APPDIR.
cat > "$appdir/AppRun" <<'EOF'
#!/bin/sh
APPDIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
exec "$APPDIR/OITOOLS" "$@"
EOF
chmod +x "$appdir/AppRun"

cat > "$appdir/oitools.desktop" <<'EOF'
[Desktop Entry]
Type=Application
Name=OITOOLS
GenericName=Optical Interferometry
Comment=Explore, model and reconstruct optical interferometry data
Exec=AppRun %f
Icon=oitools
Terminal=false
Categories=Science;Astronomy;Education;
MimeType=application/x-oifits;
StartupWMClass=OITOOLS
EOF

# install-desktop.sh belongs to the tarball, not here: an AppImage is registered by the
# desktop itself (or appimaged), and a script writing absolute paths into ~/.local would point
# at a mount that disappears when the application exits.
rm -f "$appdir/install-desktop.sh"

# ── build ────────────────────────────────────────────────────────────────────
# Each option reaches mksquashfs as its own --mksquashfs-opt; the tool does not split a
# quoted string.
opts=(--mksquashfs-opt -b --mksquashfs-opt "$block")
case $comp in
    zstd|gzip) opts+=(--mksquashfs-opt -Xcompression-level --mksquashfs-opt "$level") ;;
esac

echo "packing $(du -sh "$appdir" | cut -f1) with $comp, ${block} blocks, level $level ..."
ARCH=x86_64 "$tool" --appimage-extract-and-run --comp "$comp" "${opts[@]}" "$appdir" "$out"
chmod +x "$out"
printf 'done: %s  (%s)\n' "$out" "$(du -h "$out" | cut -f1)"
