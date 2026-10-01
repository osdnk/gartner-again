"""Validate multi-move identities and reproduce their mathematical comparisons."""

from decimal import Decimal as D, localcontext
from functools import lru_cache
from math import ceil
from random import Random

from verify import functions, masses, optimal_rate
from check_protocol import matrix_product, multiply


@lru_cache(maxsize=None)
def line(alpha, offset):
    weights = masses(alpha, offset, ceil(35 * alpha + abs(offset) + 12))
    return functions(weights, 0, D(1))


def elementary(values):
    coefficients = [D(1)]
    for value in values:
        coefficients.append(D(0))
        for j in range(len(coefficients) - 1, 0, -1):
            coefficients[j] += value * coefficients[j - 1]
    return coefficients


def orthogonal_probabilities(alpha, point, repetition):
    lines = [line(alpha, coordinate) for coordinate in point]
    ratios = [entry[2] for entry in lines]
    even = sum(elementary(ratios)[::2])
    result = []
    for i, (minus, plus, _) in enumerate(lines):
        coefficients = elementary(ratios[:i] + ratios[i + 1:])
        weight = sum(coefficients[j] / (j + 1) for j in range(0, len(coefficients), 2))
        result.append((minus * weight / (repetition * even), plus * weight / (repetition * even)))
    return result


def rotate(vector, exponent, degree):
    result = []
    for start in range(0, len(vector), degree):
        block = [0] * degree
        for j, value in enumerate(vector[start:start + degree]):
            target = j + exponent
            block[target % degree] = value * (-1 if (target // degree) % 2 else 1)
        result.extend(block)
    return result


def inner(left, right):
    return sum(x * y for x, y in zip(left, right))


def unit_probabilities(alpha, point, directions, repetition):
    norm2 = inner(directions[0], directions[0])
    lines = [line(alpha, D(inner(point, v)) / norm2) for v in directions]
    count = len(directions)
    probabilities = []
    for i, (minus, plus, ratio) in enumerate(lines):
        other_ratio = lines[(i + count // 2) % count][2]
        scale = D(2) / (count * repetition * (1 + ratio * other_ratio))
        probabilities.append((minus * scale, plus * scale))
    return probabilities


def fold(poly, directions):
    block = len(poly) // directions
    return [sum(poly[j::block]) % 2 for j in range(block)]


def main():
    with localcontext() as context:
        context.prec = 100
        tolerance = D("1e-50")
        rng = Random(20260910)
        worst_orthogonal = D(0)
        worst_unit = D(0)
        orthogonal_count = 0
        unit_count = 0
        for alpha in map(D, ("0.4", "0.8", "1.2", "2")):
            two = optimal_rate(alpha)
            delta = (two - 1) / (two + 1)
            for count in (1, 2, 4, 8):
                optimum = (1 + delta**count) / (1 - delta**count)
                for _ in range(12):
                    point = tuple(D(rng.randrange(-12, 13)) / 8 for _ in range(count))
                    for repetition in (optimum, 2 * optimum):
                        probabilities = orthogonal_probabilities(alpha, point, repetition)
                        assert min(p for pair in probabilities for p in pair) >= -tolerance
                        assert sum(sum(pair) for pair in probabilities) <= 1 + tolerance
                        incoming = D(0)
                        for i in range(count):
                            for sign, column in ((-1, 0), (1, 1)):
                                source = list(point)
                                source[i] -= sign
                                probability = orthogonal_probabilities(alpha, tuple(source), repetition)[i][column]
                                density_ratio = (-(source[i]**2 - point[i]**2) / (2 * alpha**2)).exp()
                                incoming += density_ratio * probability
                        residual = abs(incoming * repetition - 1)
                        assert residual < tolerance, (alpha, count, point, residual)
                        worst_orthogonal = max(worst_orthogonal, residual)
                        orthogonal_count += 1

            for count in (2, 4, 8):
                degree = 16
                for _ in range(5):
                    vector = [rng.randrange(-2, 3) for _ in range(2 * degree)]
                    vector[0] = 1
                    norm2 = inner(vector, vector)
                    directions = [rotate(vector, j * degree // count, degree) for j in range(count)]
                    assert all(inner(v, v) == norm2 for v in directions)
                    assert all(inner(directions[j], directions[j + count // 2]) == 0 for j in range(count // 2))
                    point = [rng.randrange(-2, 3) for _ in vector]
                    repetition = (two + 1 / two) / 2
                    probabilities = unit_probabilities(alpha, point, directions, repetition)
                    assert min(p for pair in probabilities for p in pair) >= -tolerance
                    assert sum(sum(pair) for pair in probabilities) <= 1 + tolerance
                    incoming = D(0)
                    for i, direction in enumerate(directions):
                        for sign, column in ((-1, 0), (1, 1)):
                            source = [x - sign * v for x, v in zip(point, direction)]
                            probability = unit_probabilities(alpha, source, directions, repetition)[i][column]
                            density_ratio = (-D(inner(source, source) - inner(point, point)) / (2 * alpha**2 * norm2)).exp()
                            incoming += density_ratio * probability
                    residual = abs(incoming * repetition - 1)
                    assert residual < tolerance, (alpha, count, residual)
                    worst_unit = max(worst_unit, residual)
                    unit_count += 1

        ring_count = 0
        for degree in (16, 32):
            for count in (2, 4, 8):
                for _ in range(8):
                    ell, m, modulus = 2, 2, 12289
                    f = rng.choice((1, 3, 5))
                    secret = [[f] + [0] * (degree - 1)] + [[rng.randrange(-2, 3) for _ in range(degree)] for _ in range(ell + m)]
                    a0 = [[[rng.randrange(modulus) for _ in range(degree)] for _ in range(ell)] for _ in range(m)]
                    a0s = matrix_product(a0, secret[1:ell + 1], modulus)
                    public_b = [[pow(f, -1, modulus) * (x + e) % modulus for x, e in zip(row, error)] for row, error in zip(a0s, secret[ell + 1:])]
                    matrix = []
                    for i in range(m):
                        first = [-2 * x for x in public_b[i]]
                        first[0] += modulus * int(i == 0)
                        identity = [[2 * int(i == j)] + [0] * (degree - 1) for j in range(m)]
                        matrix.append([first] + [[2 * x for x in p] for p in a0[i]] + identity)
                    block = degree // count
                    challenge = [rng.randrange(2) for _ in range(block)] + [0] * (degree - block)
                    lifted = [0] * degree
                    for j in range(block):
                        if challenge[j]:
                            lifted[j + block * rng.randrange(count)] = rng.choice((-1, 1))
                    assert fold(lifted, count) == challenge[:block]
                    y = [[rng.randrange(-30, 31) for _ in range(degree)] for _ in secret]
                    z = [[x + shift for x, shift in zip(row, multiply(lifted, s))] for row, s in zip(y, secret)]
                    assert matrix_product(matrix, y, modulus) == matrix_product(matrix, z, modulus)
                    assert [(x - c) % 2 for x, c in zip(fold(z[0], count), challenge)] == fold(y[0], count)
                    ring_count += 1

        print(f"Orthogonal identities: {orthogonal_count} point/factor tests, worst relative residual {worst_orthogonal:.3E}")
        print(f"Unit-family identities: {unit_count} tests, worst relative residual {worst_unit:.3E}")
        print(f"Folded signature equations: {ring_count} ring transcripts")
        print("alpha, rejection bits for 2/4/8/16 mutually orthogonal moves")
        for alpha in map(D, ("0.8", "1", "1.2", "1.5")):
            two = optimal_rate(alpha)
            delta = (two - 1) / (two + 1)
            bits = [-(2 * delta**a / (1 + delta**a)).ln() / D(2).ln() for a in (1, 2, 4, 8)]
            print(alpha, *(f"{value:.5f}" for value in bits))
        print("N, multiplier C sufficient for norm-abort bound 2^-192")
        for dimension in (1024, 1536, 2048):
            lo, hi = D(1), D(2)
            for _ in range(100):
                mid = (lo + hi) / 2
                exponent = D(dimension) * (mid**2 - 1 - 2 * mid.ln()) / (2 * D(2).ln())
                if exponent >= 192:
                    hi = mid
                else:
                    lo = mid
            print(dimension, f"{hi:.9f}")
        for dimension, multiplier in ((1024, "1.379769"), (1536, "1.307468"), (2048, "1.264879")):
            multiplier = D(multiplier)
            exponent = D(dimension) * (multiplier**2 - 1 - 2 * multiplier.ln()) / (2 * D(2).ln())
            assert exponent >= 192


if __name__ == "__main__":
    main()
