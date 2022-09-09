#!/usr/bin/env bash

set -euxo pipefail

main() {
  git fetch --prune origin "+refs/tags/*:refs/tags/*"
  bump2version --current-version $(git describe) \
    --tag --tag-name {new_version} --tag-message 'Bump version to {new_version}' patch
  git push --tags

  exit 0
}

main "$@"
