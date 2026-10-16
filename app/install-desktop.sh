#!/bin/sh
#
# Register this bundle with the desktop: menu entry, icon, and .oifits association.
#
# Run once after unpacking, from anywhere:   ./install-desktop.sh
#
# The paths are computed from where THIS script sits, so the bundle can be unpacked anywhere
# and moved afterwards by re-running it. Everything is written under ~/.local/share, so no
# administrator rights are needed and an uninstall is three `rm`s (printed at the end).
set -eu
here=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)
apps="$HOME/.local/share/applications"
icons="$HOME/.local/share/icons/hicolor/256x256/apps"
mime="$HOME/.local/share/mime/packages"
mkdir -p "$apps" "$icons" "$mime"

# share/oitools mirrors the repository's own relative paths -- that is what makes
# `OITOOLS.resource` a plain root substitution -- so the icon keeps its `app/assets` prefix.
cp "$here/share/oitools/app/assets/oitools.png" "$icons/oitools.png"

# A MIME type of our own: OIFITS is FITS underneath, so without this the desktop hands .oifits
# to whatever claims application/fits, and "Open with" never offers this application.
cat > "$mime/oitools.xml" <<EOF
<?xml version="1.0" encoding="UTF-8"?>
<mime-info xmlns="http://www.freedesktop.org/standards/shared-mime-info">
  <mime-type type="application/x-oifits">
    <comment>Optical interferometry data (OIFITS)</comment>
    <sub-class-of type="application/fits"/>
    <glob pattern="*.oifits"/>
    <icon name="oitools"/>
  </mime-type>
</mime-info>
EOF

# Terminal=false matters: with it true the desktop spawns a terminal that serves no purpose,
# since stdout goes nowhere useful and failures are written to the crash log instead.
cat > "$apps/oitools.desktop" <<EOF
[Desktop Entry]
Type=Application
Name=OITOOLS
GenericName=Optical Interferometry
Comment=Explore, model and reconstruct optical interferometry data
Exec=$here/OITOOLS %f
Icon=oitools
Terminal=false
Categories=Science;Astronomy;Education;
MimeType=application/x-oifits;
StartupWMClass=OITOOLS
EOF

command -v update-desktop-database >/dev/null 2>&1 && update-desktop-database "$apps" || true
command -v update-mime-database    >/dev/null 2>&1 && update-mime-database "$HOME/.local/share/mime" || true
command -v gtk-update-icon-cache   >/dev/null 2>&1 && gtk-update-icon-cache -f -t "$HOME/.local/share/icons/hicolor" 2>/dev/null || true

echo "Installed. OITOOLS is in the application menu, and .oifits files open with it."
echo "To undo:  rm $apps/oitools.desktop $icons/oitools.png $mime/oitools.xml"
