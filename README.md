# susyewkxsections

Holds useful ROOT scripts to fit/read 13 TeV ATLAS, CMS EWK direct gaugino production cross-sections.

##Currently works for:

* Chargino pair-production                    : C1C1 (grid)
* Chargino-neutralino associated production   : C1N2 (grid) [can be split into C1pN2 and C1mN2]
* Neutralino-neutralino associated production : N1N2 (grid) - only in hino
* Pure wino gauginos     : wino                    (comp)
* Pure higgsino gauginos : hino                    (comp)

## Usage

If the desired ROOT file grid_comp_13TeV.root doesn't exist, first run the script to do the fits and store the outputs as:

root -l 'fit_gaugino.C(grid,comp)' etc.

then to read the cross-sections simply do:

root -l 'get_gaugino.C(grid,comp,mass)' where mass is the sparticle mass in GeV

In case of any problems, please contact the ATLAS/CMS SUSY EWK conveners, happy hunting! :)

## References

These cross-sections are tabulated at:

[LHC SUSY Cross-sections Twiki (under Electroweak Sector)](https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections#Cross_sections_for_various_S_AN2)
