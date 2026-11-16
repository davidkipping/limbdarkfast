"""Print the hybrid laws' constants for one version, from python/limbdarkfast.py.

    python tools/law_constants.py cpp|fortran|julia|idl

For each law: the number of weights n_w, the number of poles K, the poles
eps (padded to 3), the matrix G (n_w x (3 + K), padded to 5 x 6) giving
each shape as a combination of the generators [mu^0, mu^2, mu^4,
(mu^2 + eps_k)^-2...], and norms (padded to 5), the unocculted flux of
each shape. Values are written with 17 significant digits, which
round-trip to the same doubles.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "python"))
import limbdarkfast  # noqa: E402

LAWS = ("hybrid2", "hybrid4", "hybrid5")


def padded(name):
    L = limbdarkfast.get_law(name)
    K = len(L.eps)
    eps = np.zeros(3)
    eps[:K] = L.eps
    G = np.zeros((5, 6))
    G[:L.n_w, :3 + K] = L.G
    norms = np.zeros(5)
    norms[:L.n_w] = L.norms
    return L.n_w, K, eps, G, norms


def num(x, lang):
    s = f"{x:.17g}"
    if lang == "fortran":
        s = (s.replace("e", "d") if "e" in s else s + "d0")
    elif lang == "idl":
        s = (s.replace("e", "d") if "e" in s else s + "d")
    elif "." not in s and "e" not in s:
        s += ".0"
    return s


def main(lang):
    for name in LAWS:
        n_w, K, eps, G, norms = padded(name)
        if lang == "cpp":
            rows = ",\n      ".join(
                "{" + ", ".join(num(v, lang) for v in row) + "}" for row in G)
            print(f"  static const LawData {name} = {{\n    {n_w}, {K},\n"
                  f"    {{{', '.join(num(v, lang) for v in eps)}}},\n"
                  f"    {{{rows}}},\n"
                  f"    {{{', '.join(num(v, lang) for v in norms)}}}}};")
        elif lang == "fortran":
            def wrap(vals, per=3):
                parts = [num(v, lang) for v in vals]
                lines = [", ".join(parts[i:i + per]) for i in range(0, len(parts), per)]
                return "[ &\n        " + ", &\n        ".join(lines) + "]"
            print(f"  ! {name}: G row by row (shape j over generators k)")
            print(f"  type(law_t), parameter :: law_{name} = law_t({n_w}, {K}, &\n"
                  f"      {wrap(eps)}, &\n"
                  f"      reshape({wrap(G.ravel())}, [5, 6], order=[2, 1]), &\n"
                  f"      {wrap(norms)})")
        elif lang == "julia":
            tup = lambda vals: "(" + ", ".join(num(v, lang) for v in vals) + ")"
            rows = ",\n         ".join(tup(row) for row in G)
            print(f"const {name.upper()} = Law(:{name}, {n_w}, {K},\n"
                  f"    {tup(eps)},\n"
                  f"    ({rows}),\n"
                  f"    {tup(norms)})")
        elif lang == "idl":
            arr = lambda vals: "[" + ", ".join(num(v, lang) for v in vals) + "]"
            rows = ", $\n          ".join(", ".join(num(v, lang) for v in row) for row in G)
            print(f"    '{name}': return, {{nw: {n_w}, K: {K}, $\n"
                  f"        eps: {arr(eps)}, $\n"
                  f"        G: reform([ $\n          {rows}], 6, 5), $\n"
                  f"        norms: {arr(norms)}}}")
        else:
            print(f"# {name}: n_w={n_w} K={K}")
            print("eps   =", [num(v, lang) for v in eps])
            for row in G:
                print("G row =", [num(v, lang) for v in row])
            print("norms =", [num(v, lang) for v in norms])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "cpp")
