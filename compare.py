"""Reproduce the original parameter estimates and evaluate the improved rule."""

from decimal import Decimal as D, localcontext, ROUND_CEILING

from verify import arctan_inverse, optimal_rate


def original_rate(alpha, pi):
    q = (-pi**2 * alpha**2 / 2).exp()
    return 1 + 2 * alpha * (2 * pi).sqrt() * q / ((-1 / (2 * alpha**2)).exp() * (1 - q**4))


def width_for_rate(target, bound, old_width):
    lo, hi = old_width / 2, old_width
    for _ in range(100):
        mid = (lo + hi) / 2
        if optimal_rate(mid / bound) > target:
            lo = mid
        else:
            hi = mid
    return hi


def entropy_bytes(width, degree, ell, m, compression, pi):
    first = D(degree * (ell + 1)) * (2 * pi * D(1).exp() * width**2).ln()
    second = D(degree * m) * (8 * pi * D(1).exp() * (width / compression)**2).ln()
    return int(((first + second) / (16 * D(2).ln())).to_integral_value(rounding=ROUND_CEILING)) + 64


def norm_bounds(width, degree, ell, m, compression):
    signing = (D("1.02") * width * D(degree * (ell + m + 1)).sqrt()).to_integral_value(rounding=ROUND_CEILING)
    verifying = (signing + D(degree * m).sqrt() * (D(compression) / 4 + 1)).to_integral_value(rounding=ROUND_CEILING)
    return signing, verifying


def main():
    with localcontext() as ctx:
        ctx.prec = 90
        pi = 16 * arctan_inverse(D(5)) - 4 * arctan_inverse(D(239))
        rows = [
            ("I", 12289, 1, 2, "110.07", 128, 58, 256, 775, 4178, 5649),
            ("II", 50177, 3, 2, "48.00", 55, 80, 128, 1184, 2199, 2946),
            ("III", 50177, 4, 3, "79.60", 95, 128, 128, 1694, 4386, 5301),
        ]
        for name, modulus, ell, m, printed_bound, width, kappa, compression, expected_size, expected_bs, expected_bv in rows:
            degree = 256
            width = D(width)
            bound = D(printed_bound)
            alpha = width / bound
            original = original_rate(alpha, pi)
            improved = optimal_rate(alpha)
            new_width = width_for_rate(original, bound, width)
            assert abs(optimal_rate(new_width / bound) - original) < D("1e-28")
            old_size = entropy_bytes(width, degree, ell, m, compression, pi)
            assert old_size == expected_size
            assert norm_bounds(width, degree, ell, m, compression) == (expected_bs, expected_bv)
            listed_bv = {"I": 5649, "II": 2946, "III": 5300}[name]
            assert expected_bv - listed_bv == int(name == "III")
            listed_rejection = {"I": D("0.460"), "II": D("0.642"), "III": D("0.617")}[name]
            assert (1 - original**(-kappa)).quantize(D("0.001")) == listed_rejection
            assert 32 + (degree*m*(modulus-1).bit_length()+7)//8 == {"I": 928, "II": 1056, "III": 1568}[name]
            new_size = entropy_bytes(new_width, degree, ell, m, compression, pi)
            print(f"Set {name}: Bk={bound:.10f}, alpha={alpha:.10f}")
            for label, factor in (("original", original), ("improved", improved)):
                rejection = 1 - 1 / factor
                full_rejection = 1 - factor ** (-kappa)
                print(f"  {label}: M={factor:.12f}, step=2^({rejection.ln()/D(2).ln():.6f}), stage=2^({full_rejection.ln()/D(2).ln():.6f}), attempts={factor**kappa:.9f}")
            print(f"  attempt reduction percent={100*(1-(improved/original)**kappa):.4f}")
            print(f"  fixed-rate width={width} -> {new_width:.9f}, reduction={100*(1-new_width/width):.4f}%")
            print(f"  entropy-size bytes={old_size} -> {new_size}, saved={old_size-new_size}, percentage={100*D(old_size-new_size)/old_size:.4f}")
            print(f"  norm bounds={norm_bounds(new_width,degree,ell,m,compression)}")
            print(f"  public key bytes={32 + (degree*m*(modulus-1).bit_length()+7)//8}")


if __name__ == "__main__":
    main()
