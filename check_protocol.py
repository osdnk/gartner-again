"""Check rounding and protocol identities, not cryptographic security."""

from random import Random


def centred(value, modulus):
    return (value + modulus // 2) % modulus - modulus // 2


def high_bits(value, modulus, tau):
    return tau * (((value % (2 * modulus)) + tau // 2) // tau) % (2 * (modulus - 1))


def multiply(left, right):
    degree = len(left)
    result = [0] * degree
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            if i + j < degree:
                result[i + j] += x * y
            else:
                result[i + j - degree] -= x * y
    return result


def matrix_product(matrix, vector, modulus):
    return [[sum(parts) % modulus for parts in zip(*(multiply(a, b) for a, b in zip(row, vector)))] for row in matrix]


def main():
    scalar_count = 0
    for modulus, tau in ((12289, 256), (50177, 128)):
        images = set()
        for high in range(0, 2 * (modulus - 1), tau):
            for bit in (0, 1):
                image = (high + bit) % (2 * modulus)
                assert image not in images
                images.add(image)
        for w in range(2 * modulus):
            high = high_bits(w, modulus, tau)
            low = centred(w - high, 2 * modulus)
            bit = w % 2
            assert -tau // 2 - 2 <= low <= tau // 2 - 1
            assert low % 2 == bit
            assert abs((bit - low) // 2) <= tau // 4 + 1
            for response in (-modulus - 3, -2, 0, 7, modulus + 4):
                partial = (w - 2 * response) % (2 * modulus)
                hint = (high - high_bits(partial, modulus, tau)) % (2 * (modulus - 1))
                reconstructed_high = (high_bits(partial, modulus, tau) + hint) % (2 * (modulus - 1))
                numerator = reconstructed_high - partial + bit
                assert numerator % 2 == 0
                reconstructed = centred(numerator // 2, modulus)
                assert reconstructed_high == high
                assert reconstructed == centred(response + (bit - low) // 2, modulus)
                assert abs(reconstructed) <= abs(response) + tau // 4 + 1
            scalar_count += 1

    rng = Random(20260909)
    transcript_count = 0
    for degree in (8, 16):
        for ell, m, modulus, tau in ((1, 2, 12289, 256), (3, 2, 50177, 128), (4, 3, 50177, 128)):
            for _ in range(8):
                # Odd constant f gives an invertible special case of the key equation.
                f = rng.choice((1, 3, 5))
                secret = [[f] + [0] * (degree - 1)] + [[rng.randrange(-3, 4) for _ in range(degree)] for _ in range(ell + m)]
                a0 = [[[rng.randrange(modulus) for _ in range(degree)] for _ in range(ell)] for _ in range(m)]
                a0s = matrix_product(a0, secret[1:ell + 1], modulus)
                inverse_f = pow(f, -1, modulus)
                b = [[inverse_f * (x + e) % modulus for x, e in zip(row, error)] for row, error in zip(a0s, secret[ell + 1:])]
                matrix = []
                for i in range(m):
                    first = [-2 * x for x in b[i]]
                    first[0] += modulus if i == 0 else 0
                    identity = [[2 * int(i == j)] + [0] * (degree - 1) for j in range(m)]
                    matrix.append([first] + [[2 * x for x in p] for p in a0[i]] + identity)
                image = matrix_product(matrix, secret, 2 * modulus)
                assert image == [[modulus * int(i == 0)] + [0] * (degree - 1) for i in range(m)]
                challenge = [rng.randrange(2) for _ in range(degree)]
                signed_challenge = [bit * rng.choice((-1, 1)) for bit in challenge]
                y = [[rng.randrange(-100, 101) for _ in range(degree)] for _ in secret]
                z = [[x + shift for x, shift in zip(row, multiply(signed_challenge, s))] for row, s in zip(y, secret)]
                w = matrix_product(matrix, y, 2 * modulus)
                az = matrix_product(matrix, z, 2 * modulus)
                a1z1 = matrix_product([row[:ell + 1] for row in matrix], z[:ell + 1], 2 * modulus)
                for i in range(m):
                    for j in range(degree):
                        qc = modulus * challenge[j] if i == 0 else 0
                        assert (az[i][j] - qc) % (2 * modulus) == w[i][j]
                        partial = (a1z1[i][j] - qc) % (2 * modulus)
                        assert partial == (w[i][j] - 2 * z[ell + 1 + i][j]) % (2 * modulus)
                        bit = (z[0][j] - challenge[j]) % 2 if i == 0 else 0
                        assert bit == w[i][j] % 2
                        high = high_bits(w[i][j], modulus, tau)
                        reconstructed = centred((high - partial + bit) // 2, modulus)
                        assert (a1z1[i][j] + 2 * reconstructed - qc) % (2 * modulus) == (high + bit) % (2 * modulus)
                transcript_count += 1
    print(f"Checked {scalar_count} scalar residues and {transcript_count} ring transcripts.")


if __name__ == "__main__":
    main()
