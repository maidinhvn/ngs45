"""ngs45 — recover the 45S nrDNA transcribed unit from Illumina short reads.

Short reads cannot span a full rDNA repeat, so ngs45 *assembles* the unit:
it baits rDNA reads against a conserved 45S seed, assembles them with SPAdes,
resolves a single repeat unit (monomer) out of the tandem-repeat graph, orients
and annotates it, then maps reads back to quantify ribotype heterogeneity.
"""

from __future__ import annotations

__version__ = "0.1.0"
