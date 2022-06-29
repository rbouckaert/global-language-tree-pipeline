## Monos package

This is a [BEAST](http://beast2.org) package to generate timing constraints based on posterior distributions of independent analyses represented by posterior tree samples.


## How to build

Clone the project.

Monos uses the following packages:

* beast2
* BEASTLabs
* Babel
and whatever packages are required to run the input XML file (see below).

To build the package, run `ant`. This builds a file in `build/dist/monos.addon.v1.0.0.zip`

## How to install

To install by hand: go to BEAST package directory (e.g. `~/.beast/2.6` on Linux, `C:\Users\<yourname>\BEAST` on Windows).

Create a folder `monos`

Change to the `monos` directory and unzip `build/dist/monos.addon.v1.0.0.zip`.

Reset the class path for BEAST (start BEAUti, select menu `File/Clear class path`.

## How to generate XMLs


Use the `applauncher` utility that comes with BEAST with the `GenerateLexicalConstraints` argument. For examples from within the monos data directory:

'applauncher GenerateLexicalConstraints -glottologTree dplace/iso.tree -xml dplace/geo-rc1197+almostnewwals+corrcals4.xml -treeID Tree.t:DPLACE -treeConfig treeset.cfg -burnin 10 -out /tmp/dplace.xml'

For dplace, the data (including configuration files) can be downloaded [here](https://github.com/rbouckaert/monos/releases/download/v0.0.1/monos-data.tgz).


`GenerateLexicalConstraints` has the following options:

* glottologTree 	Newick file containing the monophyletic constraints from Glottolog
* xml 	file name of BEAST XML file containing the model for which to create constraints for
* treeID [string]	name of the tree in the XML file that contains all taxa. This taxon set should contain all taxa in the target output.
* treeConfig [filename]	configuration file that identifies lexical analyses. This should be a tab delimited file with three columns: name of the posterior tree set, name of the taxon map file, and scale factor for trees used as multiplier for the tree height of trees in the tree set.
* burnin [integer]	percentage of trees to used as burn-in (and will be ignored)
* out 	output file. Print to stdout if not specified
* threshold [double]	threshold above which clade support in lexical analyses is deemed high enough to justify adding a monophyly constraint in the pruned tree.

## How to run XMLs

The XMLs require BEAST 2.6 to be installed with the AlmostDistributions and GEO_SPHERE packages. 
The GEO_SPHERE package is in the default package repository it is easiest to start BEAUti (a program that is part of BEAST), and select the menu File/Manage packages. A package manager dialog pops up where you can select the GEO_SPHERE package and click the Install/Upgrade button. The BEASTLabs package should be automatically installed as well.

If the AlmostDistributions package is listed as well, just click on it to select it, and hit the Install/Upgrade button. If the AlmostDistributions package is not listed, you may need to add a package repository by clicking the "Package repositories" button. A window pops up where you can click "Add URL" and add "https://raw.githubusercontent.com/CompEvol/CBAN/master/packages-extra.xml" in the entry. The AlmostDistributions package should now be listed in the package manager dialog.

After that, running BEAST on the XML files should work. Be aware that the working directory should contain the `kml` folder containing the two KML files with geographical constraints.



