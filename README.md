# Global language tree pipeline

This repository contains the code and supporting data for the preprint:

Bouckaert, R., Redding, D., Trebski, A., Sheehan, O., Kyritsis, T., Gray, R., Jones, K. E., & Atkinson, Q. (2026). *Global language diversification is linked to socio-ecology and threat status*. SocArXiv. [https://doi.org/10.31235/osf.io/f8tr6_v2](https://doi.org/10.31235/osf.io/f8tr6_v2)

Please cite the paper when making use of the code or derived data in this repository.

## Overview

The project has two main components:

- `monos/`: a BEAST 2 package used to generate XML files and constraints for the global language phylogeny workflow.
- `TreeSetAnalysisScripts/`: the downstream R analysis workflow used on the posterior treeset, including socio-ecological regressions, threat analyses, EDGE calculations, and manuscript figures/tables.

## Where to start

If you want to work with the downstream manuscript analyses, start here:

- [TreeSetAnalysisScripts/README.md](TreeSetAnalysisScripts/README.md): main TreeSet workflow README, including system requirements, installation, demo, usage, and reproduction notes.
- [TreeSetAnalysisScripts/demo/README.md](TreeSetAnalysisScripts/demo/README.md): focused instructions for the lightweight three-step demo workflow.

If you want to work with the BEAST-side tree generation tools, start here:

- [monos/README.md](monos/README.md)

## Data availability

Large inputs are not fully included in a standard `git clone`.

- GitHub Releases for this repository contain large tree and analysis assets used by the TreeSet workflow:
[https://github.com/rbouckaert/global-language-tree-pipeline/releases](https://github.com/rbouckaert/global-language-tree-pipeline/releases)
- Additional large environmental/spatial inputs and exact unpacking locations are documented in [TreeSetAnalysisScripts/README.md](TreeSetAnalysisScripts/README.md).

## Citation

If you use this repository, please cite the preprint above and, where relevant, the specific software/component documentation in:

- [TreeSetAnalysisScripts/README.md](TreeSetAnalysisScripts/README.md)
- [monos/README.md](monos/README.md)

