#!/bin/sh
# Sync the public snapshot (analysis code + frozen ledgers) to the public repository.
# The private development repo stays authoritative; this mirrors exactly the paths below.
set -e
SRC="$(cd "$(dirname "$0")" && pwd)"
DST="${1:-$HOME/Documents/Github/entropic-superfluid-gravity}"
mkdir -p "$DST"
EXTRA=""
[ -f "$SRC/.sync_exclude" ] && EXTRA="--exclude-from=$SRC/.sync_exclude"
rsync -a --delete --exclude '__pycache__' $EXTRA "$SRC/analysis" "$SRC/derived" "$DST/"
for f in README.md LICENSE CITATION.cff sync_public.sh; do cp "$SRC/$f" "$DST/$f"; done
printf '__pycache__/\n.DS_Store\n' > "$DST/.gitignore"
echo "synced -> $DST (commit and push from there)"
