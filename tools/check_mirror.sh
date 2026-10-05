#!/bin/sh
# Compare OmicOne's R/ against scStudio's, ignoring the package-name rename.
#
# The two repos share their single-cell code; only the multi-omics shell and the
# WES pipeline are meant to differ. A shared file that shows up here is a sync
# that was missed -- which is exactly how a call to `scstudio_theme()` survived
# in OmicOne's fct_overview.R and broke the Import overview at runtime.
#
# Usage:  tools/check_mirror.sh [path-to-scStudio]
# Exit:   0 = only the expected files differ, 1 = an unexpected divergence.

set -eu

OMIC=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SC=${1:-"$OMIC/../scStudio"}

if [ ! -d "$SC/R" ]; then
  echo "scStudio not found at $SC" >&2
  exit 2
fi

# Files that are allowed to differ, or to exist in one repo only (plus every
# mod_wes_*.R / fct_wes*.R: the WES pipeline is OmicOne's own).
EXPECTED="app_landing.R app_server.R app_ui.R steps.R mod_placeholder.R
mod_report.R fct_wes.R mod_wes_clin.R mod_wes_compare.R mod_wes_driver.R
mod_wes_hetero.R mod_wes_import.R mod_wes_lolli.R mod_wes_onco.R mod_wes_sig.R
mod_wes_summary.R mod_wes_surv.R mod_wes_titv.R mod_wes_tmb.R"

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/o" "$tmp/s"

normalise() {
  perl -pe 's/OmicOne/XSTUDIOX/g; s/omicone/xstudiox/g;
            s/scStudio/XSTUDIOX/g;  s/scstudio/xstudiox/g' "$1"
}
for f in "$OMIC"/R/*.R; do normalise "$f" > "$tmp/o/$(basename "$f")"; done
for f in "$SC"/R/*.R;   do normalise "$f" > "$tmp/s/$(basename "$f")"; done

status=0
diff -rq "$tmp/o" "$tmp/s" 2>/dev/null | while read -r line; do
  file=$(printf '%s' "$line" | sed 's/.*[ /]\([A-Za-z0-9_.]*\.R\).*/\1/')
  case "$file" in mod_wes_*|fct_wes*) continue ;; esac   # the WES pipeline is OmicOne's own
  case " $(echo $EXPECTED) " in
    *" $file "*) ;;
    *) echo "UNEXPECTED DIVERGENCE: $file" ;;
  esac
done > "$tmp/report"

if [ -s "$tmp/report" ]; then
  cat "$tmp/report"
  echo
  echo "Shared files must stay identical apart from the package name."
  echo "Sync with:  perl -pe 's/OmicOne/scStudio/g; s/omicone/scstudio/g'"
  status=1
else
  echo "Mirror OK: only the expected files differ."
fi
exit $status
