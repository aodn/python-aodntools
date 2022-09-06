#!/usr/bin/env bash

set -euxo pipefail

main() {
  git fetch --prune origin "+refs/tags/*:refs/tags/*"
  NEW_VERSION=$(bump2version --current-version $(git describe) --list --no-commit --tag \
  --tag-message 'Bump version to {new_version}' patch | grep -oP '^new_version=\K.*$')
  git push origin tag v$NEW_VERSION

  exit 0
}

main "$@"
