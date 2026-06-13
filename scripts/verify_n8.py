# -*- coding: utf-8 -*-
"""
Numerical verification (n=8, r=3) for
"Partner-Midword Joint Weight Enumerators of NP1CCs and Diamond Codes".

Checks:
  Weight gap: Type I diff=1, Type II diff in {0,2}
  Rigidity:  M_k = 2 Q_{k-1} + R_{k-1} + R_{k+1}
  Marginal:  A_w = P1_w + P1_{w-1} + Q_w + Q_{w-2} + 2 R_w
  Support:   Type A <=> Q=R=M=0 ;  Type B <=> P1=0
  Diamond:   A_i = D_i + 2 D_{i-1} + D_{i-2};  D = W/(1+t)^2
  Type A closed form:  P1_i = ((n-i) B_i + (i+1) B_{i+1})/n         (B = ext. perfect code wd)
  Type B closed form:  Q_i = ((n-1-i)/(n-1)) B'_i ,  R_i = (i/(n-1)) B'_i   (B' = perfect code wd)
Plus a small Type-C search over coordinate permutations / translates.
"""

from itertools import combinations

def wt(x):           return bin(x).count('1')
def dist(a, b):      return wt(a ^ b)

# ---------- length-7 Hamming perfect code ----------
def syndrome7(x):
    s = 0
    for j in range(7):
        if (x >> j) & 1:
            s ^= (j + 1)        # column of coordinate j is (j+1) in 3 bits
    return s

H7 = [x for x in range(1 << 7) if syndrome7(x) == 0]     # 16 codewords, min dist 3
assert len(H7) == 16

def extend(x, pos):              # append overall parity at bit position `pos`
    return x | ((wt(x) & 1) << pos)

eH8 = [extend(x, 7) for x in H7]                          # extended Hamming [8,16,4]

# ---------- partner / midword analysis ----------
def analyze(code, n):
    code = set(code)
    cw = sorted(code)
    # partner = unique codeword within distance 2
    partner, pdist = {}, {}
    ok_partner = True
    for c in cw:
        near = [d for d in cw if d != c and dist(c, d) <= 2]
        if len(near) != 1:
            ok_partner = False
            partner[c] = None
        else:
            partner[c] = near[0]
            pdist[c] = dist(c, near[0])
    # symmetry check
    for c in cw:
        if partner[c] is not None and partner.get(partner[c]) != c:
            ok_partner = False

    P1, Q, R, M = {}, {}, {}, {}
    typeI = typeII = 0
    gap_ok = True
    seen = set()
    for c in cw:
        d = partner[c]
        if d is None or c in seen:
            continue
        seen.add(c); seen.add(d)
        wc, wd = sorted((wt(c), wt(d)))
        dd = dist(c, d)
        if dd == 1:
            typeI += 1
            if wd - wc != 1: gap_ok = False
            P1[wc] = P1.get(wc, 0) + 1
        elif dd == 2:
            typeII += 1
            if wd - wc not in (0, 2): gap_ok = False
            diff = c ^ d
            bits = [j for j in range(n) if (diff >> j) & 1]      # the two differing coords
            mids = [c ^ (1 << bits[0]), c ^ (1 << bits[1])]
            for m in mids:
                M[wt(m)] = M.get(wt(m), 0) + 1
                if m in code:                                    # midwords must be non-codewords
                    gap_ok = False
            if wd - wc == 2:
                Q[wc] = Q.get(wc, 0) + 1
            else:
                R[wc] = R.get(wc, 0) + 1

    A = {}
    for c in cw:
        A[wt(c)] = A.get(wt(c), 0) + 1

    # type
    has_I  = typeI  > 0
    has_II = typeII > 0
    typ = 'A' if (has_I and not has_II) else 'B' if (has_II and not has_I) else 'C'
    return dict(ok_partner=ok_partner, gap_ok=gap_ok, type=typ,
                P1=P1, Q=Q, R=R, M=M, A=A, nI=typeI, nII=typeII)

def vec(d, n):       # dict -> list 0..n
    return [d.get(i, 0) for i in range(n + 1)]

def check_propA(r, n):
    for k in range(n + 1):
        lhs = r['M'].get(k, 0)
        rhs = 2 * r['Q'].get(k - 1, 0) + r['R'].get(k - 1, 0) + r['R'].get(k + 1, 0)
        if lhs != rhs:
            return False
    return True

def check_propB(r, n):
    for w in range(n + 1):
        lhs = r['A'].get(w, 0)
        rhs = (r['P1'].get(w, 0) + r['P1'].get(w - 1, 0)
               + r['Q'].get(w, 0) + r['Q'].get(w - 2, 0) + 2 * r['R'].get(w, 0))
        if lhs != rhs:
            return False
    return True

# ---------- build the two reference codes ----------
e1 = 1 << 0
C_A = set(eH8) | set(c ^ e1 for c in eH8)                          # Type A, n=8
C_B = set(H7) | set((c ^ e1) | (1 << 7) for c in H7)              # Type B, n=8

print("=" * 64)
print("REFERENCE CODES (n = 8)")
print("=" * 64)
for name, code in [("C_A", C_A), ("C_B", C_B)]:
    r = analyze(code, 8)
    print(f"\n[{name}]  |C|={len(code)}  type={r['type']}  "
          f"(#TypeI={r['nI']}, #TypeII={r['nII']})")
    print(f"   well-defined partners : {r['ok_partner']}")
    print(f"   weight gap             : {r['gap_ok']}")
    print(f"   rigidity identity      : {check_propA(r, 8)}")
    print(f"   marginal identity      : {check_propB(r, 8)}")
    print(f"   A_w   = {vec(r['A'], 8)}")
    print(f"   P1_i  = {vec(r['P1'], 8)}")
    print(f"   Q_i   = {vec(r['Q'], 8)}")
    print(f"   R_i   = {vec(r['R'], 8)}")
    print(f"   M_k   = {vec(r['M'], 8)}")

# ---------- support separation ----------
rA, rB = analyze(C_A, 8), analyze(C_B, 8)
propD_A = (sum(rA['Q'].values()) == 0 and sum(rA['R'].values()) == 0
           and sum(rA['M'].values()) == 0 and rA['type'] == 'A')
propD_B = (sum(rB['P1'].values()) == 0 and rB['type'] == 'B')
print("\n" + "=" * 64)
print("SUPPORT SEPARATION")
print("=" * 64)
print(f"   Type A  <=>  Q=R=M=0      : {propD_A}")
print(f"   Type B  <=>  P1=0         : {propD_B}")

# ---------- Type A / Type B closed forms ----------
def wd(code, n):
    d = {}
    for c in code:
        d[wt(c)] = d.get(wt(c), 0) + 1
    return d

B  = wd(eH8, 8)          # extended Hamming weight distribution
Bp = wd(H7, 7)           # Hamming weight distribution
n  = 8

P1_formula = [((n - i) * B.get(i, 0) + (i + 1) * B.get(i + 1, 0)) // n for i in range(n + 1)]
Q_formula  = [((n - 1 - i) * Bp.get(i, 0)) // (n - 1) for i in range(n + 1)]
R_formula  = [(i * Bp.get(i, 0)) // (n - 1) for i in range(n + 1)]

print("\n" + "=" * 64)
print("CLOSED-FORM CROSS-CHECKS")
print("=" * 64)
print(f"   ext.Hamming  B_i  = {vec(B, 8)}")
print(f"   Hamming      B'_i = {vec(Bp, 7)}")
print(f"   Type A:  P1 formula == brute : "
      f"{P1_formula == vec(rA['P1'], 8)}   {P1_formula}")
print(f"   Type B:  Q  formula == brute : "
      f"{Q_formula == vec(rB['Q'], 8)}   {Q_formula}")
print(f"   Type B:  R  formula == brute : "
      f"{R_formula == vec(rB['R'], 8)}   {R_formula}")

# ---------- Theorem C : diamond code ----------
def build_diamond(np1cc, n):
    """Extend NP1CC (length n) to ENP1CC (length n+1) then add midwords."""
    Cstar = [extend(c, n) for c in np1cc]                 # length n+1, all even weight
    code = set(Cstar)
    mids = set()
    seen = set()
    info = analyze(Cstar, n + 1)
    if not info['ok_partner']:
        return None
    cwset = set(Cstar)
    part = {}
    for c in Cstar:
        near = [d for d in cwset if d != c and dist(c, d) <= 2]
        part[c] = near[0]
    for c in Cstar:
        d = part[c]
        if c in seen: continue
        seen.add(c); seen.add(d)
        if dist(c, d) != 2:
            return None
        bits = [j for j in range(n + 1) if ((c ^ d) >> j) & 1]
        mids.add(c ^ (1 << bits[0]))
        mids.add(c ^ (1 << bits[1]))
    return set(Cstar) | mids

def poly_div_1pt2(A):
    """divide sum A_i t^i by (1+t)^2 ; return D or None if not divisible over Z."""
    # divide by (1+t) twice
    def div1(p):
        q = [0] * len(p)
        rem = list(p)
        for i in range(len(p) - 1):
            q[i] = rem[i]
            rem[i + 1] -= q[i]
            rem[i] = 0
        if rem[-1] != 0:
            return None
        return q[:-1] if len(q) > 1 else q
    p1 = div1(A)
    if p1 is None: return None
    p2 = div1(p1)
    return p2

print("\n" + "=" * 64)
print("Theorem C  (diamond code from C_B)")
print("=" * 64)
dia = build_diamond(C_B, 8)
if dia is None:
    print("   diamond construction FAILED (not a clean diamond)")
else:
    rd = analyze(dia, 9)
    Ad = vec(rd['A'], 9)
    # each codeword must have exactly two neighbours at distance 1 (quotient matrix s11=2)
    deg1 = all(sum(1 for y in dia if y != x and dist(x, y) == 1) == 2 for x in dia)
    D = poly_div_1pt2(Ad)
    # check A_i = D_i + 2 D_{i-1} + D_{i-2}
    recon_ok = D is not None and all(
        Ad[i] == (D[i] if i < len(D) else 0)
                 + 2 * (D[i - 1] if 0 <= i - 1 < len(D) else 0)
                 + (D[i - 2] if 0 <= i - 2 < len(D) else 0)
        for i in range(len(Ad)))
    print(f"   |diamond|={len(dia)}  (expect 64)")
    print(f"   every codeword has exactly 2 dist-1 neighbours : {deg1}")
    print(f"   A_i(diamond)         = {Ad}")
    print(f"   D_i = W/(1+t)^2      = {D}")
    print(f"   D non-negative & integer                       : "
          f"{D is not None and all(d >= 0 for d in D)}")
    print(f"   A_i = D_i+2D_(i-1)+D_(i-2) reconstruction      : {recon_ok}")

# ---------- small Type-C search ----------
print("\n" + "=" * 64)
print("TYPE-C SEARCH (coordinate permutations & translates, n=8)")
print("=" * 64)
import itertools, random

def permute(code, perm, n):
    out = set()
    for c in code:
        y = 0
        for j in range(n):
            if (c >> j) & 1:
                y |= 1 << perm[j]
        out.add(y)
    return out

types = {'A': 0, 'B': 0, 'C': 0, 'bad': 0}
spectra = {}        # (type, A_tuple) -> set of (P1,Q,R) tuples
random.seed(0)
trials = []
# family: translates of C_A and C_B, plus random coordinate permutations
base = [("C_A", C_A), ("C_B", C_B)]
for _ in range(40):
    perm = list(range(8)); random.shuffle(perm)
    base.append((f"permA", permute(C_A, perm, 8)))
    perm = list(range(8)); random.shuffle(perm)
    base.append((f"permB", permute(C_B, perm, 8)))
# also translates of C_B by every weight-<=2 vector
for t in [0, 1, 3, (1 | (1 << 7)), (1 << 7)]:
    base.append((f"transB", set(c ^ t for c in C_B)))

for name, code in base:
    if len(code) != 32:
        types['bad'] += 1; continue
    r = analyze(code, 8)
    if not r['ok_partner']:
        types['bad'] += 1; continue
    types[r['type']] += 1
    key = (r['type'], tuple(vec(r['A'], 8)))
    spec = (tuple(vec(r['P1'], 8)), tuple(vec(r['Q'], 8)), tuple(vec(r['R'], 8)))
    spectra.setdefault(key, set()).add(spec)

print(f"   codes scanned: {len(base)}   type counts: "
      f"A={types['A']} B={types['B']} C={types['C']} (invalid={types['bad']})")
# strict-refinement witness: same (type,A) but different joint spectrum
witness = [(k, v) for k, v in spectra.items() if len(v) > 1]
if witness:
    print("   STRICT-REFINEMENT WITNESS FOUND (same A, different joint spectrum):")
    for k, v in witness:
        print(f"      type={k[0]} A={list(k[1])}  -> {len(v)} distinct spectra")
else:
    print("   no strict-refinement witness among these families")
    print("   (Type C absent at n=8 under linear constructions; length-7 perfect")
    print("    codes are all Hamming-equivalent. The length-16 script supplies")
    print("    the Type C strict-refinement examples.)")

print("\nDONE.")
