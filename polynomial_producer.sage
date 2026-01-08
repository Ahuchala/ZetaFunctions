# Pairs and their indices: 0↦01, 1↦02, ..., 9↦34
pairs = [(0,1), (0,2), (0,3), (0,4),
         (1,2), (1,3), (1,4),
         (2,3), (2,4),
         (3,4)]

pair_index = {p: i for i, p in enumerate(pairs)}

def induced_on_pairs(cycle_str):
    """
    Input:  cycle_str – a permutation of S_5 in cycle notation,
            acting on {1,2,3,4,5}.
            (Your 0,1,2,3,4 correspond to 1,2,3,4,5 here.)

    Output: list [σ̄(0), ..., σ̄(9)] giving the induced permutation
            on the 10 pairs in one-line notation (0..9).
    """
    S5 = SymmetricGroup(5)
    sigma1 = S5(cycle_str)        # S_5 element on {1,...,5}

    # Convert to a map on {0,...,4} (shift down by 1)
    def sigma(x):
        return int(sigma1(x+1)) - 1

    images = []
    for (a, b) in pairs:
        a1, b1 = sigma(a), sigma(b)
        if a1 > b1:
            a1, b1 = b1, a1
        images.append(pair_index[(a1, b1)])
    return images

induced_on_pairs('(1,2,3)')
