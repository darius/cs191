"""
Let's play with chemical reaction networks.
"""

import itertools

p2 = """\
x + y -> 2z
z + x -> 2x
z + y -> 2y
""".splitlines()

## parse_CRN(p2)
#. x + y -> 2z
#. x + z -> 2x
#. y + z -> 2y

## print parse_CRN(p2).rate_equation()
#. dx/dt = -1*(x y) + (x z)
#. dy/dt = -1*(x y) + (y z)
#. dz/dt = 2*(x y) + -1*(x z) + -1*(y z)
#. 

example = """\
    a <-> 2b
a + c <-> d
    d  -> b + e
b + e  -> a + c
""".splitlines()

## parse_CRN(example)
#. a -> 2b
#. a + c -> d
#. b + e -> a + c
#. 2b -> a
#. d -> a + c
#. d -> b + e
## parse_CRN(example).species()
#. set(['a', 'c', 'b', 'e', 'd'])
## print parse_CRN(example).rate_equation()
#. da/dt = -1*(a c) + -1*(a) + (b e) + (b^2) + (d)
#. db/dt = 2*(a) + -1*(b e) + -2*(b^2) + (d)
#. dc/dt = -1*(a c) + (b e) + (d)
#. dd/dt = (a c) + -2*(d)
#. de/dt = -1*(b e) + (d)
#. 

"""
Rate equation for example:

da/dt =    -k0 a - k1 a c + k2 b e +   k3 b^2 + k4 d
db/dt =   2 k0 a          - k2 b e - 2 k3 b^2        + k5 d

So a src complex including (sym, n) adds -n*k*(factor) to the RHS for
d(sym)/dt, where factor includes sym^n.
In other words:

   d(sym)/dt = ... + -k*n*(...sym^n...) + ...

A dst complex including (sym, n) adds n*k*(factor) to the RHS for
d(sym)/dt.
"""

class CRN:
    def __init__(self, reactions):
        self.reactions = sorted(reactions)
    def __repr__(self):
        return '\n'.join(map(show_reaction, self.reactions))
    def species(self):
        return set(symbol
                   for src, dst in self.reactions
                   for symbol, n in src + dst)
    def rate_equation(self):
        rhses = {}
        def add(sym, factor, coeff):
            terms = rhses.setdefault(sym, {})
            terms[factor] = terms.get(factor, 0) + coeff
        for src, dst in self.reactions:
            k = 1  # rate constant for this reaction -- XXX fill in
            factor = '(%s)' % (' '.join(pow(sym, n)
                                        for sym, n in sorted(src)))
            for sym, n in src: add(sym, factor, -n * k)
            for sym, n in dst: add(sym, factor, n * k)
        return '\n'.join('d%s/dt = %s' % (sym, total(rhs))
                         for sym, rhs in sorted(rhses.items()))

def total(terms): 
    return ' + '.join(mul(factor, coeff)
                      for factor, coeff in sorted(terms.items())
                      if coeff != 0)
def mul(s, n): return s if n == 1 else '%s*%s' % (n, s)
def pow(s, n): return s if n == 1 else '%s^%s' % (s, n)

def show_reaction((src, dst)):
    return '%s -> %s' % (show_complex(src), show_complex(dst))

def show_complex(complex):
    return ' + '.join('%s%s' % ('' if n == 1 else n, species)
                      for species, n in complex)

def parse_CRN(equations):
    return CRN(reaction
               for equation in equations
               for reaction in parse_equation(equation))

def parse_equation(s):
    arrow = '<->' if '<->' in s else '->'
    complexes = chop(s, arrow)
    src, dst = [parse_complex(complex) for complex in complexes]
    if arrow == '<->':
        return [(src, dst), (dst, src)]
    else:
        return [(src, dst)]

def chop(s, sep):
    assert sep in s, "%r not in %r" % (sep, s)
    return s.split(sep, 1)

def parse_complex(complex):
    return tuple(sorted([parse_term(term) for term in complex.split('+')]))

def parse_term(s):
    s = s.strip()
    coeff   = ''.join(itertools.takewhile(lambda c: c.isdigit(), s))
    species = ''.join(itertools.dropwhile(lambda c: c.isdigit(), s))
    assert species
    return species, (int(coeff) if coeff else 1)
