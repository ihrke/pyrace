echo "Starting to update \n"

#go to home and setup git
cd $HOME
git config --global user.email "travis@travis-ci.org"
git config --global user.name "Travis"

#using token clone gh-pages branch
git clone --quiet --branch=notebook_tests https://${GH_TOKEN}@github.com/ihrke/pyrace.git  notebook_tests > /dev/null

#go into diractory and copy data we're interested in to that directory
cd notebook_tests
python rerun_nb.py

#add, commit and push files
git commit -a -m "Travis build $TRAVIS_BUILD_NUMBER"
git push -fq origin notebook_tests > /dev/null
echo 'pushed'

