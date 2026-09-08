#!/bin/bash
#
# format.sh -- apply the repository code style (.clang-format) to C sources.
#
# Usage:
#   ./format.sh                 format lines changed vs. the base branch
#   ./format.sh <ref>           format lines changed vs. <ref>
#   ./format.sh --staged        format lines currently staged for commit
#   ./format.sh --check         report unformatted changed lines, change nothing
#   ./format.sh --all           format every tracked C source file
#
# The default mode only touches lines you have edited, so the rest of the tree
# and its git-blame history are left alone.

set -euo pipefail

BASE_BRANCH=${FORMAT_BASE:-origin/main_fixes}
MIN_VERSION=18

cd "$(git rev-parse --show-toplevel)"

# -- locate a new enough clang-format ----------------------------------------

find_binary()
{
  local candidate
  for candidate in clang-format-19 clang-format-18 clang-format; do
    if command -v "$candidate" >/dev/null 2>&1; then
      echo "$candidate"
      return 0
    fi
  done
  return 1
}

CF=$(find_binary) || {
  echo "error: clang-format not found. Install clang-format $MIN_VERSION or newer:" >&2
  echo "         Debian/Ubuntu: sudo apt install clang-format-$MIN_VERSION" >&2
  echo "         macOS:         brew install clang-format" >&2
  exit 1
}

CF_VERSION=$("$CF" --version | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -1)
CF_MAJOR=${CF_VERSION%%.*}

if [ "$CF_MAJOR" -lt "$MIN_VERSION" ]; then
  echo "error: found $CF version $CF_VERSION, but this repository needs $MIN_VERSION or newer." >&2
  echo "       Older versions produce different output and would fight with other contributors." >&2
  exit 1
fi

# clang-format-diff ships next to clang-format, sometimes only version-suffixed.
CFD=""
for candidate in "clang-format-diff-$CF_MAJOR" clang-format-diff; do
  if command -v "$candidate" >/dev/null 2>&1; then
    CFD=$candidate
    break
  fi
done

# -- file selection -----------------------------------------------------------

# Vendored and third-party code we do not restyle. The profile_util submodule is
# excluded automatically because git ls-files does not descend into submodules.
EXCLUDE_RE='^src/extern/'

tracked_sources()
{
  git ls-files '*.c' '*.h' | grep -Ev "$EXCLUDE_RE"
}

# -- modes --------------------------------------------------------------------

format_diff()   # $1 = git diff arguments
{
  if [ -z "$CFD" ]; then
    echo "error: clang-format-diff not found; incremental formatting needs it." >&2
    echo "       Use './format.sh --all' or install the clang-format-$CF_MAJOR package." >&2
    exit 1
  fi

  # -U0 so only edited lines are considered; -p1 strips the a/ b/ diff prefixes.
  git diff -U0 --no-color "$@" -- '*.c' '*.h' \
    | grep -v "^+++ b/src/extern/" \
    | "$CFD" -p1 -i -style=file -binary "$CF"
}

check_diff()
{
  local out
  out=$(git diff -U0 --no-color "$@" -- '*.c' '*.h' | "$CFD" -p1 -style=file -binary "$CF")
  if [ -n "$out" ]; then
    echo "$out"
    echo >&2
    echo "error: the lines above are not formatted. Run ./format.sh to fix them." >&2
    exit 1
  fi
  echo "All changed lines are correctly formatted."
}

case "${1:---default}" in
  --all)
    echo "Formatting all tracked C sources..."
    tracked_sources | xargs "$CF" -i -style=file
    echo "Done. Review with 'git diff' before committing."
    ;;
  --staged)
    format_diff --cached
    ;;
  --check)
    check_diff "$BASE_BRANCH"
    ;;
  --help | -h)
    sed -n '2,15p' "$0" | sed 's/^# \?//'
    ;;
  --default)
    git rev-parse --verify --quiet "$BASE_BRANCH" >/dev/null \
      || { echo "error: base ref '$BASE_BRANCH' not found. Pass a ref, or set FORMAT_BASE." >&2; exit 1; }
    format_diff "$BASE_BRANCH"
    ;;
  *)
    git rev-parse --verify --quiet "$1" >/dev/null \
      || { echo "error: '$1' is not a git ref. See ./format.sh --help." >&2; exit 1; }
    format_diff "$1"
    ;;
esac