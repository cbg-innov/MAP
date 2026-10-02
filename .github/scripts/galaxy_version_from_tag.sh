#!/usr/bin/env bash
# On a release tag, set the Galaxy wrapper's version from the tag, so a MAP release
# needs no edit to macros.xml:
#     vX.Y.Z           -> TOOL_VERSION X.Y.Z, VERSION_SUFFIX 0  (new MAP release)
#     vX.Y.Z+galaxyN   -> TOOL_VERSION X.Y.Z, VERSION_SUFFIX N  (wrapper-only change)
# TOOL_VERSION also selects the container (ghcr.io/cbg-innov/map:X.Y.Z), which must
# already be on GHCR. On branch pushes and pull requests this does nothing, and the
# values committed in macros.xml are tested as they are.
set -euo pipefail

macros=galaxy/map/macros.xml

if [[ "${GITHUB_REF:-}" != refs/tags/v* ]]; then
    echo "Not a release tag: using macros.xml as committed."
    exit 0
fi

tag="${GITHUB_REF_NAME#v}"
if [[ ! "$tag" =~ ^([0-9]+\.[0-9]+\.[0-9]+)(\+galaxy([0-9]+))?$ ]]; then
    echo "::error::Tag $GITHUB_REF_NAME is not vX.Y.Z or vX.Y.Z+galaxyN."
    exit 1
fi
version="${BASH_REMATCH[1]}"
suffix="${BASH_REMATCH[3]:-0}"

image="ghcr.io/cbg-innov/map:$version"
if ! docker manifest inspect "$image" > /dev/null 2>&1; then
    echo "::error::$image is not on GHCR. Push the MAP image before pushing tag $GITHUB_REF_NAME."
    exit 1
fi

sed -i \
    -e "s|<token name=\"@TOOL_VERSION@\">[^<]*</token>|<token name=\"@TOOL_VERSION@\">$version</token>|" \
    -e "s|<token name=\"@VERSION_SUFFIX@\">[^<]*</token>|<token name=\"@VERSION_SUFFIX@\">$suffix</token>|" \
    "$macros"
grep -q "<token name=\"@TOOL_VERSION@\">$version</token>" "$macros"
grep -q "<token name=\"@VERSION_SUFFIX@\">$suffix</token>" "$macros"
echo "Galaxy wrapper version set to $version+galaxy$suffix (container $image)."
