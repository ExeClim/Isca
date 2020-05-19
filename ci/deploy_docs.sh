#!/bin/bash
# exit with nonzero exit code if anything fails
set -e

# Clone *this* git repo, but only the gh-pages branch.
echo Cloning gh-pages...
if [[ ! -d $GH_PAGES_DIR ]]; then
    git clone -q -b gh-pages --single-branch https://${GITHUB_API_KEY}@github.com/${TRAVIS_REPO_SLUG}.git $GH_PAGES_DIR
fi
cd $GH_PAGES_DIR

# inside this git repo we pretend to be isca
git config user.name "isca-model"
git config user.email "iscagcm@gmail.com"

# if [[ "${TRAVIS_TAG}" != "" ]]; then
#     # export VERSION=${TRAVIS_TAG%.*}
#     export VERSION=${TRAVIS_TAG}
# else
export VERSION="dev"
# fi

# The first and only commit to this new Git repo contains all the
# files present with the commit message "Deploy to GitHub Pages".
echo Updating $VERSION docs...
rm -rf ${VERSION}
mkdir -p ${VERSION}
cp -R ${TRAVIS_BUILD_DIR}/docs/_build/html ${VERSION}/
cp -R ${TRAVIS_BUILD_DIR}/docs/_build/doctrees ${VERSION}/
touch .nojekyll
# if [[ "${VERSION}" != "dev" ]]; then
ln -snf ${VERSION} latest
# fi

# Generate our json list of versions
# echo Generating versions.json...
# ${TRAVIS_BUILD_DIR}/ci/gen_versions_json.py

echo Staging...
git add -A .
# if [[ "${VERSION}" == "dev" ]]; then
##  && `git log -1 --format='%s'` == *"dev"* ]]; then
git commit --amend --reset-author --no-edit
# else
#     git commit -m "Deploy ${VERSION} to GitHub Pages"
# fi

# Push up to gh-pages
echo Pushing...
git push --force -q https://${GITHUB_API_KEY}@github.com/${TRAVIS_REPO_SLUG}.git gh-pages
