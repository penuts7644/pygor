#!/bin/bash

LAST_VERSION="$(git rev-list --tags --skip=1 --max-count=1)"
VERSION="$(git describe --abbrev=0 --tags)"

# Generate CHANGELOG.md
changelog -m  \
  $(cut -d "/" -f1 <<< $TRAVIS_REPO_SLUG) \
  $(cut -d "/" -f2 <<< $TRAVIS_REPO_SLUG)

body="$(cat CHANGELOG.md)"

# Overwrite CHANGELOG.md with JSON data for GitHub API
jq -n \
  --arg body "$body" \
  --arg name "$VERSION" \
  --arg tag_name "$VERSION" \
  --arg target_commitish "$GIT_BRANCH" \
  '{
    body: $body,
    name: $name,
    tag_name: $tag_name,
    target_commitish: $target_commitish,
    draft: false,
    prerelease: false
  }' > CHANGELOG.md

echo "Create release $VERSION for repo: $TRAVIS_REPO_SLUG, branch: $GIT_BRANCH"
curl -H "Authorization: token $GITHUB_TOKEN" --data @CHANGELOG.md "https://api.github.com/repos/$TRAVIS_REPO_SLUG/releases"
