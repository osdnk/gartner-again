# Optimal Acceptance for Iterative Rejection Sampling

The paper uses the supplied IACR Communications in Cryptology template. Build with `make` and open `main.pdf`. The bibliography uses the existing cryptobib checkout in `../optimal-gartner/cryptobib`.

The paper has two contributions, explicit acceptance functions in the original three-condition framework and a matching optimality theorem for every positive Gaussian width. The original research notes are preserved in their respective directories.

The main body contains the introduction and technical overview, the sampling framework, optimal two-move acceptance, the basic signature scheme, four-move signatures, and the parameter and impact tables. Bit-dropping, full simulation and security proofs, Gaussian product calculations, and supporting parameter accounting follow the references as appendices. The introduction cites Lithium (ePrint 2026/1790) using the implementation results reported in its abstract, without claiming measured improvements to Lithium.

It includes the complete mathematical basic and compressed signature protocols, their correctness proofs, accepting-transcript simulators, and a conditional ordinary EUF-CMA theorem in the classical random-oracle model. The theorem explicitly includes the norm-conditioned key distribution. The before-and-after tables separate fixed-width rejection improvements from fixed-rate entropy-based signature-size estimates. They do not claim new concrete security levels or measured encoding sizes.

The four-move section includes optimal acceptance probabilities, basic and compressed protocols, correctness, accepting-transcript simulators, and a conditional security theorem for the folded self-target relation. A short remark explains why mixing more orthogonal pairs gives the same rejection factor while reducing challenge capacity.

The four-move comparison covers three lattice parameter sets at degrees 512 and 1024. Their challenge rules are fixed HW 58 on 256 positions, full binary on 256 positions with maximum HW 256, and fixed HW 82 on 512 positions. At two expected sampler attempts, the original-rule and four-move entropy estimates are 1892 to 1730, 2557 to 2362, and 2991 to 2732 bytes. Keys and challenge distributions are fixed within each comparison. Full parameters, module and coefficient SIS dimensions, widths, norm bounds, and unstructured lattice proxies appear in the paper.

The original degree-256 parameters remain a separate baseline with the exact 775, 1184, and 1694 byte counts and reduced weights 58, 80, and at most 128. All nine rejecting-mode BKZ/Core-SVP pairs reproduce exactly. The paper preserves the cited third verification bound of 5300 in that reproduction and separately explains why enforcing the displayed ceiling formula requires 5301.

For the new sets, an assumed classical hash-query budget of 2^64 gives challenge-search contributions 2^-129.553, 2^-192, and 2^-256.578. These are not complete forgery bounds, and the unstructured SIS costs do not establish the folded self-target assumption. No concrete encoder or constant-time signer is implemented.

Run `make verify` to check the formulas numerically and reproduce the comparison tables. The validation uses Python's standard library. It checks both move probabilities, their sum, the incoming Gaussian mass, and the centred-offset bound at local, global, and larger repetition factors. It also checks compression rounding exhaustively at both concrete moduli and verifies the protocol equations on small-ring test transcripts. These are checks of mathematical formulas, not a cryptographic sampler implementation.

The optional lattice-proxy validation uses SageMath and the lattice estimator. Run `make security ESTIMATOR=/path/to/lattice-estimator`. The tested estimator revision is `3e48ef421ec256afddb3e7d2249a77eab6e9ba12`. The script first checks the exact rejecting-mode baseline, then computes SIS at both verification bounds and the NTWE proxy for every new row.
