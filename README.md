# isimip-sensitivity-analysis
Using the results of the LakeEnsemblR calibration for the ISIMIP 3b simulations to discover more about lake model behaviour for a large set of lakes.

How do I contribute new code back to the `isimip-sensitivity-analysis` project?
==========================================================

In order to contribute to this code, we recommend the following workflow:

1.  "fork" this repository to your own personal github account

2.  clone the github repository to your computer:

    $git clone <git@github.com:{username}/isimip-sensitivity-analysis.git>

3.  modify code or add new functionality, save the code

4.  add the repository main to a remote main called "upstream"

    $cd isimip-sensitivity-analysis

    $git remote add upstream <git@github.com:aemon-j/isimip-sensitivity-analysis.git>

5.  before pushing your changes to your repository, pull in the current version of the aemon-j main:

    $git fetch upstream

6.  merge these differences with your own "main" version:

    $git merge upstream/main

7.  push your changes to your github repository, in addition to changes made by pulling in the aemon-j main:

    $git push

8.  submit a pull request to aemon-j main using your account at github.com
