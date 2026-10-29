#!/usr/bin/env bash
# On a release tag, set the Galaxy wrapper's version from the tag, so a MAP release
# needs no edit to macros.xml:
#     vX.Y.Z  -> TOOL_VERSION X.Y.Z (the wrapper's version)
# TOOL_VERSION also selects the container (ghcr.io/cbg-innov/map:X.Y.Z), which must
# already be on GHCR. On branch pushes and pull requests this does nothing, and the
# values committed in macros.xml are tested as they are.
set -euo pipefail

macros=galaxy/map/macros.xml

if [[ "${GITHUB_REF:-}" != refs/tags/v* ]]; then
    echo "Not a release tag: using macros.xml as committed."
    exit 0
fi

version="${GITHUB_REF_NAME#v}"
if [[ ! "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "::error::Tag $GITHUB_REF_NAME is not vX.Y.Z."
    exit 1
fi

image="ghcr.io/cbg-innov/map:$version"
if ! docker manifest inspect "$image" > /dev/null 2>&1; then
    echo "::error::$image is not on GHCR. Push the MAP image before pushing tag $GITHUB_REF_NAME."
    exit 1
fi

sed -i \
    -e "s|<token name=\"@TOOL_VERSION@\">[^<]*</token>|<token name=\"@TOOL_VERSION@\">$version</token>|" \
    "$macros"
grep -q "<token name=\"@TOOL_VERSION@\">$version</token>" "$macros"
echo "Galaxy wrapper version set to $version (container $image)."
