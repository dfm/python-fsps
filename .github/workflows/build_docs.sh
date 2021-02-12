#!/bin/sh -l
set -e

# Build the docs
cd docs
make dirhtml

# Update the gh-pages branch
git clone --branch=gh-pages https://github.com/$GITHUB_REPOSITORY _output
cd _output
rm -rf latest
mkdir -p latest
cp -r ../_build/dirhtml/* latest/

# Push the results to GitHub
git add latest
if git -c user.name='gh-actions' -c user.email='gh-actions' commit -m "Updated docs [ci skip]"; then
    git push --force https://x-access-token:$GITHUB_TOKEN@github.com/$GITHUB_REPOSITORY gh-pages
else
    echo "No changes"
fi
