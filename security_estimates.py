"""Reproduce the cited rejecting parameters before estimating new lattice proxies."""

import argparse
from decimal import Decimal as D, localcontext
import json
from pathlib import Path
import sys

from compare_four import PARAMETERS, evaluate
from verify import arctan_inverse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--estimator", required=True, type=Path)
    parser.add_argument("--baseline-only", action="store_true")
    args = parser.parse_args()
    sys.path.insert(0, str(args.estimator.resolve()))
    from estimator import LWE, ND, SIS
    from sage.all import floor, log, oo

    def summarise(costs):
        if not costs:
            raise RuntimeError("The estimator returned no attack costs.")
        best = min(costs.items(), key=lambda item: item[1]["rop"])
        return dict(attack=best[0], beta=int(best[1]["beta"]),
                    bits=int(floor(log(best[1]["rop"], 2))))

    def sis(rows, cols, modulus, bound):
        p = SIS.Parameters(n=rows, q=modulus, m=cols, length_bound=bound, norm=2)
        return summarise(SIS.estimate.rough(p, quiet=True, catch_exceptions=False))

    def ntwe(dimension, modulus, sigma):
        p = LWE.Parameters(n=dimension, q=modulus, Xs=ND.DiscreteGaussian(sigma),
                           Xe=ND.DiscreteGaussian(sigma), m=oo)
        return summarise(LWE.estimate.rough(p, quiet=True, catch_exceptions=False))

    baseline = (
        ("I", 1, 2, 12289, 2.60, 5649, (420, 122), (332, 96), (414, 120)),
        ("II", 3, 2, 50177, 1.00, 2946, (629, 183), (500, 146), (618, 180)),
        ("III", 4, 3, 50177, 1.50, 5300, (899, 262), (730, 213), (881, 257)),
    )
    for name, ell, m, modulus, sigma, bv, want_sis, want_difference, want_ntwe in baseline:
        ordinary = sis(256 * m, 256 * (ell + m + 1), modulus, bv)
        difference = sis(256 * m, 256 * (ell + m + 1), modulus, 2 * bv)
        key = ntwe(256 * (ell + 1) - 1, modulus, sigma)
        result = dict(set=name, kind="Gartner Table 1", sis=ordinary,
                      sis_difference=difference, ntwe=key)
        print(json.dumps(result), flush=True)
        assert (ordinary["beta"], ordinary["bits"]) == want_sis
        assert (difference["beta"], difference["bits"]) == want_difference
        assert (key["beta"], key["bits"]) == want_ntwe
    if args.baseline_only:
        return
    with localcontext() as context:
        context.prec = 90
        pi = 16 * arctan_inverse(D(5)) - 4 * arctan_inverse(D(239))
        for params in PARAMETERS:
            key = ntwe((params.ell + 1) * params.degree - 1, params.modulus, float(params.sigma))
            for row in evaluate(params, pi):
                result = dict(set=params.name, kind="unstructured proxies", rule=row["rule"],
                              bv=row["bv"], ntwe=key,
                              sis=sis(params.m * params.degree, params.dimension, params.modulus, row["bv"]),
                              sis_difference=sis(params.m * params.degree, params.dimension, params.modulus, 2 * row["bv"]))
                print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
