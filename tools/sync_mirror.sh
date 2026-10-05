#!/bin/sh
# Copy OmicOne's shared single-cell code into scStudio, renaming the package.
#
# OmicOne is the source of truth for every file the two repos share. Edit the
# shared file here, run this script, then run tools/check_mirror.sh. Files that
# are meant to differ (the multi-omics shell, the WES pipeline) are listed in
# EXPECTED and never copied; edit those in each repo by hand.
#
# Usage:  tools/sync_mirror.sh [path-to-scStudio]       (dry run: add -n)
# Exit:   0 = synced (or nothing to do), 2 = scStudio not found.

set -eu

DRY=0
if [ "${1:-}" = "-n" ]; then DRY=1; shift; fi

OMIC=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
SC=${1:-"$OMIC/../scStudio"}

if [ ! -d "$SC/R" ]; then
  echo "scStudio not found at $SC" >&2
  exit 2
fi

# Keep in step with tools/check_mirror.sh.
EXPECTED="app_landing.R app_server.R app_ui.R steps.R mod_placeholder.R
mod_report.R fct_wes.R mod_wes_clin.R mod_wes_compare.R mod_wes_driver.R
mod_wes_hetero.R mod_wes_import.R mod_wes_lolli.R mod_wes_onco.R mod_wes_sig.R
mod_wes_summary.R mod_wes_surv.R mod_wes_titv.R mod_wes_tmb.R"

rename() {
  perl -pe 's/OmicOne/scStudio/g; s/omicone/scstudio/g' "$1"
}

copy_one() {   # $1 = source file, $2 = destination file
  tmp=$(mktemp)
  rename "$1" > "$tmp"
  if [ -f "$2" ] && cmp -s "$tmp" "$2"; then rm -f "$tmp"; return 0; fi
  echo "sync: ${2#$SC/}"
  if [ "$DRY" -eq 1 ]; then rm -f "$tmp"; else mv "$tmp" "$2"; fi
}

for f in "$OMIC"/R/*.R; do
  b=$(basename "$f")
  case " $(echo $EXPECTED) " in
    *" $b "*) continue ;;
  esac
  copy_one "$f" "$SC/R/$b"
done

# The stylesheet is shared verbatim (apart from the class prefix). app.js is
# not: OmicOne's carries the landing-page card animations on top.
copy_one "$OMIC/inst/app/www/custom.css" "$SC/inst/app/www/custom.css"
# The explainer-animation engine and the single-cell scenes are shared too;
# explain-wes.js is OmicOne-only.
copy_one "$OMIC/inst/app/www/explain.js" "$SC/inst/app/www/explain.js"
copy_one "$OMIC/inst/app/www/explain-sc.js" "$SC/inst/app/www/explain-sc.js"

# Shared test files (the omics-specific ones are edited per repo).
for t in helper-fake.R test-logic.R test-survival.R test-state.R \
         test-sc-compute.R test-scop.R; do
  [ -f "$OMIC/tests/testthat/$t" ] && copy_one "$OMIC/tests/testthat/$t" "$SC/tests/testthat/$t"
done

# The animation checker (renamed like the code it checks).
copy_one "$OMIC/tools/check_explain.js" "$SC/tools/check_explain.js"

# The rules every model follows, identical in both repos (no rename: they
# name both packages on purpose).
for d in AGENTS.md CLAUDE.md; do
  if [ -f "$SC/$d" ] && cmp -s "$OMIC/$d" "$SC/$d"; then continue; fi
  echo "sync: $d"
  [ "$DRY" -eq 1 ] || cp "$OMIC/$d" "$SC/$d"
done

[ "$DRY" -eq 1 ] && echo "(dry run: nothing written)"
exit 0
