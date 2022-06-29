# `EdgeXMLGenerator`

For getting information from the GlottoLog newick tree file and related data and produce a BEAST analysis in `edge.xml`.

## What `EdgeXMLGenerator` does

* select taxa
	* identify all nodes marked as language in the glottolog tree
	* remove families from `familyExclusions.csv`
	* remove languages from 'languageExclusions.csv'
* identify locations from `languages_and_dialects_geo.csv` and 'languages_and_dialects_geo2.csv' and filter out taxa without locations
* define clades
	* one for each family
	* one for regional sets (e.g., Sahul)
	* some specials for internal nodes (e.g., Arabic)
* generate GlottoLog monophyletic constraint		
* add time calibrations on families and clade MRCAs
* add time calibrations on originates of MRCAs
* process XML template (e.g., add prior on tree stuff)

To complete the XML to run with BEAST
* run monos
* run splits

## Running `EdgeXMLGenerator`

Running `EdgeXMLGenerator` requires the <a href="https://github.com/rbouckaert/monos">monos</a> package to be installed, and can be run using the BEAST 2 `applauncher` utility, e.g., using

    /path/to/beast/bin/applauncher EdgeXMLGenerator -verbose false

`EdgeXMLGenerator` has the following options:
* `glottologtree` -- file contains GlottoLog newick tree file. Default `data/tree_glottolog_newick.txt`.
* `familyExclusions` -- text file with GlottoLog codes of families that should not be included in the analysis (e.g., sign languages). Default `data/familyExclusions.csv`.
* `languageExclusions` -- text file with GlottoLog codes of languages that should not be included in the analysis. Default `data/languageExclusions.dat`.
* `cladesDir` -- specify directory with clades over families, e.g., Americas. Default `data/clades`.
* `calibrations` -- NEXUS file with time calibrations for MRCAs of clades. Default `data/calibrations.nex`.
* `originateCalibrations` -- NEXUS file with time calibrations for orignate of MRCA of clades (including singletons). Default `data/originateCalibrations.nex`.
* `glottoLocations` -- CSV file downloaded from GlottoLog with locations. Default `data/languages_and_dialects_geo.csv`.
* `customLocations` -- CSV file with locations not in GlottoLog file. Default `data/languages_and_dialects_geo2.csv`.
* `template","template of BEAST XML for EDGE analysis", new XMLFile("data/template.xml`
* `out` -- output directory, edge.xml + other outputs for sanity checking are written here. Deafult `output`.
* `verbose` -- if true, generates svg files for all families + outputs more on screen. Default `false`.


`EdgeXMLGenerator` by default gets input from the `data` directory and writes output to the `output` directory. Relevant files are described below.

## data directory

* `tree_glottolog_newick.txt` -- downloaded from GlottoLog
* `languages_and_dialects_geo.csv` -- downloaded from GlottoLog

* `familyExclusions.csv` -- list of language families that should not be included in the analysis (because they are not tree-like, or ancient)
* `languageExclusions.dat` -- list of languages that should be excluded, such as ancient languages, colonial languages, creoles, pidgins
* `languages_and_dialects_geo2.csv` -- locations not found in languages_and_dialects_geo.csv downloaded from GlottoLog

* `calibrations.nex` -- list of calibrations on MRCAs of clades
* `originateCalibrations.nex` -- list of calibrations on originate of MRCAs of clades

* `template.xml` -- BEAST XML template (should not be changed)

* `clades:
* `Americas.fam` -- list of families making up the Americas
* `Sahul.fam` -- list of families making up the Sahul
* `RestOfWorld.fam` -- list of families making up the RestOfWorld
* `specials.csv` -- list of internal nodes that have a calibration (specified in calibrations.nex)


## output directory for `EdgeXMLGenerator`

* `edge.xml` -- BEAST XML containing everything but lexical constraints. Run the `GenerateLexicalConstraints` app to insert these.

* `languageExclusion.log` -- list of languages that could not be excluded because they are not part of the analysis. If they are listed here they should probably be removed from the `/data/languageExclusions.dat` file.
* `removedDueToLackOfLocation.log` -- list of languages that were removed because no location was listed in `/data/languages_and_dialects_geo.csv` or `/data/languages_and_dialects_geo2.csv`. These languages either should be added to `/data/languages_and_dialects_geo2.csv` if a location can be identified, or in `/data/languageExclusions.dat` if they should be omitted from the analysis (because they are ancient, creoles, or otherwise not appropriate to include).
* `excludedFamilies.log` -- list of families excluded from the analyses due to being all ancient, non tree like, or otherwise inappropriate. Any family here should be added to `/data/familyExclusions.csv` if the family should be excluded, or languages + locations should be identified so the family can be included.


### SVG files
For every regional clade, an SVG file will be generated. To make
With the `-verbose` option, this is done for every family as well. 
* `Africa.svg`
* `Americas-withoutAustronesia.svg`
* `Americas.svg`
* `RestOfWorld-withoutAustronesia.svg`
* `RestOfWorld.svg`
* `Sahul-withoutAustronesia.svg`
* `Sahul.svg`

 

## Upgrading GlottoLog

Download the latest `tree_glottolog_newick.txt` and`languages_and_dialects_geo.csv` files from GlottoLog and store them in the `data` directory.

Run `EdgeXMLGenerator` and verify the following:

* the `output/languageExclusion.log` file should contain no languages (see above how to deal with exceptions).
* the `output/excludedFamilies.log` file should contain no families (see above how to deal with exceptions).
* run with `-verbose` and check there are no outliers, in particular
    * Indo European tends to have a lot of colonial languages that should be excluded, so open `Indo-European.svg` to make sure they don't.
    * The same applies to all other family SVG files -- checking `Americas.svg`, `Sahul.svg` and `RestOfWorld.svg` may be a good starting point for this check.
    * If `Africa.svg` shows any languages outside Africa, they should be added to `/data/clades/RestOfWorld.fam` and possibly on of `/data/clades/Sahul.fam` or `/data/clades/Americas.fam`.
* to identify clades, compare the taxon set in previous `edge.xml` with the newly generated `edge.xml` file.
	* save old and new taxon sets into two text files, say `old` and `new`.
	* run `sort old >x ; sort new > y; diff x y` to identify differences between taxon sets.
