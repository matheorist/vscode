# -*- coding: utf-8 -*-
"""
Exact n=16 (r=4) verification for the manuscript.

Goal: realize Type C NP1CCs via the two-perfect-code product
      C = C1 x {0}  U  C2 x {1}, with C1 = Hamming[15], C2 = Vasil'ev[15],
and test the TRICHOTOMY:
      Type A  <=>  C1 = C2
      Type B  <=>  C1 cap C2 = {empty}
      Type C  <=>  0 < |C1 cap C2| < |C1|,
together with the joint-spectrum identity  P1_i = #{ c in C1 cap C2 : wt(c)=i }.
Also search for a strict-refinement witness (same A, different (P1,Q,R)).
"""
import numpy as np
from itertools import combinations

POP = np.array([bin(i).count('1') for i in range(1 << 16)], dtype=np.uint8)
def wt(x): return bin(x).count('1')

# ---------- length-15 Hamming perfect code (linear) ----------
def syndrome15(x):
    s = 0
    for j in range(15):
        if (x >> j) & 1:
            s ^= (j + 1)
    return s
Ham15 = [x for x in range(1 << 15) if syndrome15(x) == 0]   # 2048 codewords
assert len(Ham15) == 2048

# ---------- length-7 Hamming perfect code (for Vasil'ev seed) ----------
def syndrome7(x):
    s = 0
    for j in range(7):
        if (x >> j) & 1:
            s ^= (j + 1)
    return s
Ham7 = [x for x in range(1 << 7) if syndrome7(x) == 0]       # 16 codewords
par = lambda y: bin(y).count('1') & 1

def vasilev15(lam):
    """Vasil'ev code length 15 from Ham7 and lambda: Ham7 -> F2 with lam[0]=0.
       codeword bits: [ y(0..6) | (y^c)(7..13) | par(y)^lam(c) (14) ]."""
    V = []
    for c in Ham7:
        lc = lam(c)
        for y in range(1 << 7):
            w = y | ((y ^ c) << 7) | ((par(y) ^ lc) << 14)
            V.append(w)
    return V

def is_perfect15(code):
    """min distance >= 3 and size 2048  => 1-perfect (sphere-packing equality)."""
    if len(set(code)) != 2048:
        return False
    a = np.array(sorted(set(code)), dtype=np.uint16)
    for i in range(0, len(a), 256):                     # blockwise min-distance>=3
        blk = a[i:i+256]
        d = POP[np.bitwise_xor(blk[:, None], a[None, :]).astype(np.uint16)]
        d[d == 0] = 99
        if d.min() < 3:
            return False
    return True

# ---------- generic n=16 partner / spectrum analyzer ----------
def analyze16(codewords):
    a = np.array(sorted(set(codewords)), dtype=np.uint16)
    N = len(a)
    partner = np.empty(N, dtype=np.int64)
    for i in range(N):
        d = POP[np.bitwise_xor(a[i], a).astype(np.uint16)]
        idx = np.where((d >= 1) & (d <= 2))[0]
        if idx.size != 1:
            return None
        partner[i] = idx[0]
    # symmetry
    if not all(partner[partner[i]] == i for i in range(N)):
        return None
    P1, Q, R, M, A = {}, {}, {}, {}, {}
    seen = set()
    nI = nII = 0
    for i in range(N):
        A[wt(int(a[i]))] = A.get(wt(int(a[i])), 0) + 1
    for i in range(N):
        j = partner[i]
        if i in seen:
            continue
        seen.add(i); seen.add(j)
        ci, cj = int(a[i]), int(a[j])
        wi, wj = sorted((wt(ci), wt(cj)))
        dd = POP[ci ^ cj]
        if dd == 1:
            nI += 1
            P1[wi] = P1.get(wi, 0) + 1
        else:
            nII += 1
            bits = [b for b in range(16) if ((ci ^ cj) >> b) & 1]
            for m in (ci ^ (1 << bits[0]), ci ^ (1 << bits[1])):
                M[wt(m)] = M.get(wt(m), 0) + 1
            if wj - wi == 2:
                Q[wi] = Q.get(wi, 0) + 1
            else:
                R[wi] = R.get(wi, 0) + 1
    typ = 'A' if (nI and not nII) else 'B' if (nII and not nI) else 'C'
    return dict(type=typ, P1=P1, Q=Q, R=R, M=M, A=A, nI=nI, nII=nII)

def vec(d, n=16): return tuple(d.get(i, 0) for i in range(n + 1))

def product_code(C1, C2):
    return [c for c in C1] + [(c | (1 << 15)) for c in C2]

def checkA(r):                      # marginal identity
    return all(r['A'].get(w, 0) ==
               r['P1'].get(w, 0) + r['P1'].get(w-1, 0) +
               r['Q'].get(w, 0) + r['Q'].get(w-2, 0) + 2*r['R'].get(w, 0)
               for w in range(17))
def checkRig(r):                    # rigidity identity
    return all(r['M'].get(k, 0) == 2*r['Q'].get(k-1, 0) +
               r['R'].get(k-1, 0) + r['R'].get(k+1, 0) for k in range(17))

# ---------- lambda families (lam[0]=0 always) ----------
def lam_zero(c):  return 0
def lam_b(a, b): return (lambda c: ((c >> a) & 1) & ((c >> b) & 1))   # nonlinear
def lam_b3(a, b, d): return (lambda c: ((c>>a)&1) & ((c>>b)&1) & ((c>>d)&1))

print("=" * 70)
print("n=16 TRICHOTOMY + Type C realization (C1=Hamming[15], C2=Vasil'ev[15])")
print("=" * 70)

C1 = Ham15
Hset = set(C1)
wdist = lambda S: tuple(sum(1 for c in S if wt(c) == i) for i in range(16))

# --- Type A control: C2 = C1 ---
rA = analyze16(product_code(C1, C1))
print(f"\n[A-control] C2=C1            -> type {rA['type']}  "
      f"(nI={rA['nI']}, nII={rA['nII']})  PropA={checkRig(rA)} PropB={checkA(rA)}")

# --- Type B control: C2 = C1 + v  (disjoint coset, still perfect) ---
v = 1  # weight-1 vector, not in Hamming
C2B = [c ^ v for c in C1]
rB = analyze16(product_code(C1, C2B))
interB = set(C1) & set(C2B)
print(f"[B-control] C2=C1+e1 disjoint -> type {rB['type']}  "
      f"(nI={rB['nI']}, nII={rB['nII']})  |C1capC2|={len(interB)}  "
      f"PropA={checkRig(rB)} PropB={checkA(rB)}")

# --- Type C: all products of two or three seed coordinates, plus lambda=0 ---
results = []
lam_list = [("vzero", lam_zero)]
lam_list += [(f"v({a},{b})", lam_b(a, b)) for a, b in combinations(range(7), 2)]
lam_list += [(f"v3({a},{b},{d})", lam_b3(a, b, d))
             for a, b, d in combinations(range(7), 3)]
assert len(lam_list) == 57

for name, lam in lam_list:
    V = vasilev15(lam)
    if not is_perfect15(V):
        print(f"   [{name}] NOT perfect -- skipped"); continue
    Vset = set(V)
    inter = Hset & Vset
    r = analyze16(product_code(C1, V))
    if r is None:
        print(f"   [{name}] product not a valid NP1CC (partner check failed)"); continue
    # verify P1 == weight distribution of the intersection
    p1_pred = wdist(inter)
    p1_real = vec(r['P1'])[:16]
    ok_inter = (p1_pred == p1_real)
    results.append((name, r, len(inter)))
    print(f"\n[{name}] |C1capC2|={len(inter):4d}  type={r['type']}  "
          f"(nI={r['nI']}, nII={r['nII']})")
    print(f"      PropA(rigidity)={checkRig(r)}  PropB(marginal)={checkA(r)}  "
          f"P1==wd(C1capC2): {ok_inter}")

# ---------- strict-refinement witness search ----------
print("\n" + "=" * 70)
print("STRICT-REFINEMENT WITNESS SEARCH (Type C products, same A, diff spectrum)")
print("=" * 70)
byA = {}
for name, r, ni in results:
    if r['type'] != 'C':
        continue
    key = vec(r['A'])
    spec = (vec(r['P1']), vec(r['Q']), vec(r['R']))
    byA.setdefault(key, []).append((name, spec))

assert len(byA) == 1
named_specs = {name: spec for name, spec in next(iter(byA.values()))}
assert named_specs["v(0,1)"] != named_specs["v(0,2)"]
print("  *** WITNESS: v(0,1) and v(0,2) share A but have distinct spectra ***")
print(f"      v(0,1): P1={named_specs['v(0,1)'][0]}")
print(f"      v(0,2): P1={named_specs['v(0,2)'][0]}")

# ---------- full 57-function classification ----------
classes = {}
for name, r, ni in results:
    assert r['type'] == 'C'
    assert checkRig(r) and checkA(r)
    spec = (vec(r['P1']), vec(r['Q']), vec(r['R']))
    classes.setdefault(spec, []).append((name, ni))

assert len(results) == 57
assert len(classes) == 12
assert sorted(len(members) for members in classes.values()) == \
       [1, 1, 1, 2, 2, 4, 4, 4, 6, 8, 8, 16]
assert {ni for members in classes.values() for _, ni in members} == {256, 320}

print("\n" + "=" * 70)
print("FULL 57-FUNCTION CLASSIFICATION: 12 DISTINCT JOINT SPECTRA")
print("=" * 70)
ordered = sorted(classes.items(), key=lambda item: (item[0][0][3:8], item[1][0][0]))
for index, (spec, members) in enumerate(ordered, start=1):
    p1, q, r = spec
    representative, intersection_size = members[0]
    print(f"  class {index:2d}: multiplicity={len(members):2d}, "
          f"representative={representative:12s}, |C1capC2|={intersection_size}")
    print(f"            P1[3..7]={p1[3:8]}")

print("\nDONE.")
