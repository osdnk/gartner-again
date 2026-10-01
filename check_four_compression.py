"""Check four-move rounding, compressed transcripts, and verifier rejection."""

from hashlib import shake_256
from math import ceil, isqrt, sqrt
from random import Random

from check_multimodal import fold
from check_protocol import centred, matrix_product, multiply


def high_bits(value, modulus, delta):
    return delta * (((value % modulus) + delta // 2) // delta) % (modulus - 1)


def hash_challenge(high, parity, message, degree):
    # A deterministic test oracle, not a specification of a hash-to-ring format.
    encoded = repr((high, parity, message)).encode("ascii")
    digest = shake_256(encoded).digest((degree // 2 + 7) // 8)
    bits = [(digest[j // 8] >> (j % 8)) & 1 for j in range(degree // 2)]
    return bits + [0] * (degree // 2)


def reconstruct(matrix, first, hint, modulus, delta):
    partial = matrix_product([row[:len(first)] for row in matrix], first, modulus)
    high = [[(high_bits(x, modulus, delta) + h) % (modulus - 1)
             for x, h in zip(row, hrow)] for row, hrow in zip(partial, hint)]
    last = [[centred(h - x, modulus) for h, x in zip(hrow, row)]
            for hrow, row in zip(high, partial)]
    return first + last, high


def verify(matrix, signature, message, modulus, delta, bound, degree, ell, m):
    first, hint, challenge = signature
    if len(first) != ell + 1 or len(hint) != m or len(challenge) != degree:
        return False
    if any(len(row) != degree for row in first + hint):
        return False
    if any(not isinstance(x, int) for row in first + hint for x in row):
        return False
    if any(x not in (0, 1) for x in challenge[:degree // 2]) or any(challenge[degree // 2:]):
        return False
    if any(not 0 <= x < modulus - 1 or x % delta for row in hint for x in row):
        return False
    response, high = reconstruct(matrix, first, hint, modulus, delta)
    if sum(x * x for row in response for x in row) > bound * bound:
        return False
    parity = [(x - c) % 2 for x, c in zip(fold(first[0], 2), challenge)]
    return hash_challenge(high, parity, message, degree) == challenge


def main():
    residues = 0
    for modulus, delta in ((17, 2), (17, 4), (17, 8), (17, 16),
                           (12289, 128), (50177, 64)):
        for value in range(modulus):
            high = high_bits(value, modulus, delta)
            low = centred(value - high, modulus)
            assert abs(low) <= delta // 2 + 1
            for response in (-2 * modulus - 7, -modulus, -1, 0, 1, 7, 2 * modulus + 3):
                partial = (value - response) % modulus
                hint = (high - high_bits(partial, modulus, delta)) % (modulus - 1)
                rebuilt_high = (high_bits(partial, modulus, delta) + hint) % (modulus - 1)
                rebuilt = centred(rebuilt_high - partial, modulus)
                assert rebuilt_high == high
                assert rebuilt == centred(response - low, modulus)
                assert abs(rebuilt) <= abs(response) + delta // 2 + 1
            residues += 1

    rng = Random(20260911)
    transcripts = 0
    for degree in (16, 32):
        for modulus, delta in ((12289, 128), (50177, 64)):
            for _ in range(16):
                ell, m = 3, 2
                # A nonconstant f=1+2X, with its exact inverse modulo X^d+1 and q.
                f = [1, 2] + [0] * (degree - 2)
                denominator = (1 + pow(2, degree, modulus)) % modulus
                inverse_f = [(pow(-2, j, modulus) * pow(denominator, -1, modulus)) % modulus
                             for j in range(degree)]
                assert [x % modulus for x in multiply(f, inverse_f)] == [1] + [0] * (degree - 1)
                secret = [f] + [[rng.randrange(-2, 3) for _ in range(degree)]
                                for _ in range(ell + m)]
                a0 = [[[rng.randrange(modulus) for _ in range(degree)]
                       for _ in range(ell)] for _ in range(m)]
                a0s = matrix_product(a0, secret[1:ell + 1], modulus)
                public_b = [[x % modulus for x in multiply(inverse_f, [a + e for a, e in zip(row, error)])]
                            for row, error in zip(a0s, secret[ell + 1:])]
                matrix = []
                for i in range(m):
                    identity = [[int(i == j)] + [0] * (degree - 1) for j in range(m)]
                    matrix.append([[-x for x in public_b[i]]] + a0[i] + identity)
                assert matrix_product(matrix, secret, modulus) == [[0] * degree for _ in range(m)]
                y = [[rng.randrange(-100, 101) for _ in range(degree)] for _ in secret]
                commitment = matrix_product(matrix, y, modulus)
                high = [[high_bits(x, modulus, delta) for x in row] for row in commitment]
                parity = fold(y[0], 2)
                message = f"transcript-{transcripts}"
                challenge = hash_challenge(high, parity, message, degree)
                lifted = [0] * degree
                for j, bit in enumerate(challenge[:degree // 2]):
                    if bit:
                        lifted[j + rng.randrange(2) * (degree // 2)] = rng.choice((-1, 1))
                assert fold(lifted, 2) == challenge[:degree // 2]
                assert sum(abs(x) for x in lifted) == sum(challenge)
                response = [[x + t for x, t in zip(row, multiply(lifted, s))]
                            for row, s in zip(y, secret)]
                first, last = response[:ell + 1], response[ell + 1:]
                partial = [[(w - z) % modulus for w, z in zip(wrow, zrow)]
                           for wrow, zrow in zip(commitment, last)]
                hint = [[(h - high_bits(x, modulus, delta)) % (modulus - 1)
                         for h, x in zip(hrow, row)] for hrow, row in zip(high, partial)]
                signing_bound = isqrt(sum(x * x for row in response for x in row)) + 1
                bound = ceil(signing_bound + sqrt(degree * m) * (delta // 2 + 1))
                signature = (first, hint, challenge)
                assert verify(matrix, signature, message, modulus, delta, bound, degree, ell, m)
                rebuilt, rebuilt_high = reconstruct(matrix, first, hint, modulus, delta)
                assert rebuilt_high == high
                assert matrix_product(matrix, rebuilt, modulus) == high
                assert fold(rebuilt[0], 2) == fold(response[0], 2)
                # The accepting simulator uses the same deterministic transcript map.
                assert matrix_product(matrix, response, modulus) == commitment
                assert [(x - c) % 2 for x, c in zip(fold(response[0], 2), challenge)] == parity
                bad_hint = [row[:] for row in hint]
                bad_hint[0][0] = modulus - 1
                assert not verify(matrix, (first, bad_hint, challenge), message,
                                  modulus, delta, bound, degree, ell, m)
                bad_challenge = challenge[:]
                bad_challenge[-1] = 1
                assert not verify(matrix, (first, hint, bad_challenge), message,
                                  modulus, delta, bound, degree, ell, m)
                assert not verify(matrix, signature, message, modulus, delta, 0, degree, ell, m)
                transcripts += 1
    print(f"Four-move compression: {residues} exhaustive scalar residues, {transcripts} ring transcripts.")
    print(f"Verifier rejects malformed hints, out-of-support challenges, and oversized responses in {3 * transcripts} tests.")


if __name__ == "__main__":
    main()
