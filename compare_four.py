"""Reproduce four-move comparisons across distinct lattice parameter sets."""

from dataclasses import dataclass
from decimal import Decimal as D, ROUND_CEILING, localcontext
from math import comb

from compare import entropy_bytes, norm_bounds, original_rate
from verify import arctan_inverse, optimal_rate


@dataclass(frozen=True)
class Parameters:
    name: str
    degree: int
    ell: int
    m: int
    sigma: str
    secret_bound: int
    kappa: int
    weight_mode: str
    search_target: int
    modulus: int = 50177
    delta: int = 64

    @property
    def challenge_count(self):
        positions = self.degree // 2
        return 2**positions if self.weight_mode == "full" else comb(positions, self.kappa)

    @property
    def challenge_bytes(self):
        bits = (self.challenge_count - 1).bit_length()
        return max(32, (bits + 7) // 8)

    @property
    def dimension(self):
        return self.degree * (self.ell + self.m + 1)


PARAMETERS = (
    Parameters("A", 512, 2, 2, "0.90", 58, 58, "fixed", 128),
    Parameters("B", 512, 3, 2, "1.00", 68, 256, "full", 192),
    Parameters("C", 1024, 1, 2, "1.00", 85, 82, "fixed", 256),
)

EXPECTED = {
    "A": ((1892, 3449, 4505), (1864, 3249, 4305), (1730, 2431, 3487)),
    "B": ((2557, 4918, 5974), (2529, 4676, 5732), (2362, 3460, 4516)),
    "C": ((2991, 6565, 8059), (2949, 6200, 7694), (2732, 4624, 6118)),
}


def four_rate(alpha):
    two = optimal_rate(alpha)
    return (two + 1 / two) / 2


def solve_width(rate, target, bound):
    low, high = bound / 4, 2 * bound
    assert rate(low / bound) > target > rate(high / bound)
    for _ in range(160):
        middle = (low + high) / 2
        if rate(middle / bound) > target:
            low = middle
        else:
            high = middle
    return high


def evaluate(params, pi, attempts=2):
    target = (D(attempts).ln() / params.kappa).exp()
    bound = D(params.secret_bound)
    result = []
    rates = (("original", lambda alpha: original_rate(alpha, pi)),
             ("two", optimal_rate), ("four", four_rate))
    for name, rate in rates:
        width = solve_width(rate, target, bound)
        assert abs(rate(width / bound)**params.kappa - attempts) < D("1e-40")
        size = entropy_bytes(width, params.degree, params.ell, params.m, 2 * params.delta, pi)
        size += params.challenge_bytes - 32
        bs, bv = norm_bounds(width, params.degree, params.ell, params.m, 2 * params.delta)
        result.append(dict(rule=name, width=width, size=size, bs=int(bs), bv=int(bv)))
    return result


def main():
    with localcontext() as context:
        context.prec = 90
        pi = 16 * arctan_inverse(D(5)) - 4 * arctan_inverse(D(239))
        for params in PARAMETERS:
            assert (params.modulus - 1) % params.delta == 0
            assert params.weight_mode in ("full", "fixed")
            assert 0 < params.kappa <= params.degree // 2
            assert params.challenge_count >= 2**(params.search_target + 64)
            raw_bound = D(params.sigma) * D(params.degree * (params.ell + params.m + 4)).sqrt()
            assert params.secret_bound == int(raw_bound.to_integral_value(rounding=ROUND_CEILING))
            if params.weight_mode == "fixed":
                assert comb(params.degree // 2, params.kappa - 1) < 2**(params.search_target + 64)
                assert params.challenge_count >= (2**64 + 1) * 2**params.search_target
            entropy = D(params.challenge_count).ln() / D(2).ln()
            pk_bytes = 32 + (params.degree * params.m * (params.modulus - 1).bit_length() + 7) // 8
            print(f"Set {params.name}: d={params.degree}, q={params.modulus}, ell={params.ell}, m={params.m}, sigma={params.sigma}, Bk={params.secret_bound}, Delta={params.delta}")
            print(f"  SIS module matrix={params.m} x {params.ell + params.m + 1}; coefficient matrix={params.degree * params.m} x {params.dimension}")
            print(f"  NTWE proxy dimension={(params.ell + 1) * params.degree - 1}; key bytes={pk_bytes}")
            print(f"  challenge: {params.weight_mode} HW={params.kappa}, positions={params.degree // 2}, h={entropy:.9f}, bytes={params.challenge_bytes}, query term=2^-{entropy - 64:.9f}")
            result = evaluate(params, pi)
            assert tuple((row["size"], row["bs"], row["bv"]) for row in result) == EXPECTED[params.name]
            for row in result:
                print(f"  {row['rule']}: r={row['width']:.10f}, bytes={row['size']}, Bs={row['bs']}, Bv={row['bv']}, 2Bv={2 * row['bv']}")
                assert 2 * row["bv"] < params.modulus
            assert result[2]["size"] < result[1]["size"] < result[0]["size"]
            saving = result[0]["size"] - result[2]["size"]
            print(f"  size saving: {saving} bytes ({100 * D(saving) / result[0]['size']:.6f}%)")
            fixed_width = result[0]["width"]
            for name, factor in (("original", original_rate(fixed_width / params.secret_bound, pi)),
                                 ("two", optimal_rate(fixed_width / params.secret_bound)),
                                 ("four", four_rate(fixed_width / params.secret_bound))):
                attempts = factor**params.kappa
                bits = -(1 - 1 / attempts).ln() / D(2).ln()
                print(f"  fixed-width {name}: rejection=2^-{bits:.9f}, attempts={attempts:.9f}")


if __name__ == "__main__":
    main()
