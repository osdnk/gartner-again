"""Numerical checks for the acceptance functions. Uses only the standard library."""

from decimal import Decimal as D, localcontext
from math import ceil


def arctan_inverse(n):
    x = D(1) / n
    term = x
    result = term
    k = 1
    while abs(term) > D("1e-155"):
        term *= -x * x
        result += term / (2 * k + 1)
        k += 1
    return result


def masses(alpha, offset, radius):
    return {
        k: (-(offset + k) ** 2 / (2 * alpha**2)).exp()
        for k in range(-radius, radius + 1)
    }


def functions(p, k, repetition):
    own = sum(w for j, w in p.items() if (j - k) % 2 == 0)
    other = sum(w for j, w in p.items() if (j - k) % 2 != 0)
    ratio = other / own
    right_even = sum(w for j, w in p.items() if j >= k and (j - k) % 2 == 0)
    right_odd = sum(w for j, w in p.items() if j > k and (j - k) % 2 != 0)
    left_even = sum(w for j, w in p.items() if j <= k and (j - k) % 2 == 0)
    left_odd = sum(w for j, w in p.items() if j < k and (j - k) % 2 != 0)
    f = (ratio * right_even - right_odd) / (repetition * p[k])
    g = (ratio * left_even - left_odd) / (repetition * p[k])
    return f, g, ratio


def optimal_rate(alpha):
    p = masses(alpha, D(0), ceil(40 * alpha + 10))
    even = sum(w for k, w in p.items() if k % 2 == 0)
    odd = sum(w for k, w in p.items() if k % 2 != 0)
    return even / odd


def main():
    with localcontext() as ctx:
        ctx.prec = 160
        pi = 16 * arctan_inverse(D(5)) - 4 * arctan_inverse(D(239))
        tolerance = D("1e-70")
        worst_balance = D(0)
        worst_sum = D(0)
        count = 0
        widths = list(map(D, ("0.15", "0.25", "0.4", "0.6", "0.8", "1", "1.2", "1.5", "2", "4")))
        for alpha in widths:
            global_rate = optimal_rate(alpha)
            for i in range(17):
                offset = D(i) / 16
                p = masses(alpha, offset, ceil(40 * alpha + 12))
                even = sum(w for k, w in p.items() if k % 2 == 0)
                odd = sum(w for k, w in p.items() if k % 2 != 0)
                local_rate = max(even / odd, odd / even)
                assert local_rate <= global_rate + tolerance
                # Check the optimum on each line, the common global rate,
                # and a deliberately larger admissible rate.
                for repetition in (local_rate, global_rate, 2 * global_rate):
                    radius = ceil(5 * alpha) + 1
                    values = {k: functions(p, k, repetition) for k in range(-radius - 1, radius + 2)}
                    for k in range(-radius, radius + 1):
                        f, g, ratio = values[k]
                        assert f >= -tolerance and g >= -tolerance
                        assert f + g <= 1 + tolerance
                        sum_error = abs(f + g - ratio / repetition)
                        incoming = p[k + 1] * values[k + 1][0] + p[k - 1] * values[k - 1][1]
                        balance_error = abs(incoming * repetition / p[k] - 1)
                        worst_sum = max(worst_sum, sum_error)
                        worst_balance = max(worst_balance, balance_error)
                        assert sum_error < tolerance
                        assert balance_error < tolerance
                        count += 1
        print(f"Validated {count} point/rate combinations on {len(widths) * 17} Gaussian lines.")
        print(f"Maximum relative balance residual: {worst_balance:.3E}")
        print(f"Maximum acceptance-sum residual: {worst_sum:.3E}")
        print("alpha  M_original  M_optimal  rejection_original_bits  rejection_optimal_bits  reduction_percent")
        for alpha in map(D, ("0.6", "0.8", "1", "1.2", "1.5", "2")):
            q = (-pi**2 * alpha**2 / 2).exp()
            original = 1 + 2 * alpha * (2 * pi).sqrt() * q / ((-1 / (2 * alpha**2)).exp() * (1 - q**4))
            optimal = optimal_rate(alpha)
            rej_original = 1 - 1 / original
            rej_optimal = 1 - 1 / optimal
            original_bits = -rej_original.ln() / D(2).ln()
            optimal_bits = -rej_optimal.ln() / D(2).ln()
            reduction = 100 * (1 - rej_optimal / rej_original)
            print(f"{alpha}  {original:.9f}  {optimal:.9f}  {original_bits:.5f}  {optimal_bits:.5f}  {reduction:.3f}")


if __name__ == "__main__":
    main()
