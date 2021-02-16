#!/bin/sh -l
set -e

echo "Stable: $STABLE"
VERSION=`python -c 'import fsps;print(fsps.__version__)'`
echo "Version: $VERSION"

# Build the docs
cd docs
make dirhtml

# Update the gh-pages branch
git clone --branch=gh-pages https://github.com/$GITHUB_REPOSITORY _output
cd _output
rm -rf latest
mkdir -p latest
cp -r ../_build/dirhtml/* latest/
git add latest

if [ "$STABLE" = "true" ]; then
    rm -rf $VERSION
    cp -r latest $VERSION
    git add $VERSION
fi

# Push the results to GitHub
if git -c user.name='gh-actions' -c user.email='gh-actions' commit -m "Updated docs [ci skip]"; then
    git push --force https://x-access-token:$GITHUB_TOKEN@github.com/$GITHUB_REPOSITORY gh-pages
else
    echo "No changes"
fi
