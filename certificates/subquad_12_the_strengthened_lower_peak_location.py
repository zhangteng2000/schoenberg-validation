"""Exact interval proof of the sharper lower peak location for q in [1/20,1].

All decisions use outward-rounded integer intervals.  Transcendental
enclosures use finite Taylor series and explicit rational tails.
"""
from fractions import Fraction as Q
from math import factorial
from functools import lru_cache

S = 10**70
ORDER = 60
LINEAR = Q(9, 2)
QUADRATIC = Q(23, 2)


class I:
    def __init__(self, x=0, hi=None, raw=False):
        if raw:
            self.lo, self.hi = x, hi
            return
        if isinstance(x, I):
            self.lo, self.hi = x.lo, x.hi
            return
        x = Q(x)
        self.lo = x.numerator*S//x.denominator
        self.hi = -((-x.numerator*S)//x.denominator)

    @classmethod
    def bounds(cls, a, b):
        return cls(I(a).lo, I(b).hi, True)

    def __add__(self, other):
        other = I(other)
        return I(self.lo+other.lo, self.hi+other.hi, True)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo, True)

    def __sub__(self, other):
        return self+-I(other)

    def __rsub__(self, other):
        return I(other)+-self

    def __mul__(self, other):
        other = I(other)
        vals = [self.lo*other.lo, self.lo*other.hi,
                self.hi*other.lo, self.hi*other.hi]
        return I(min(vals)//S, -((-max(vals))//S), True)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = I(other)
        assert not other.lo <= 0 <= other.hi
        vals = [Q(S*S, other.lo), Q(S*S, other.hi)]
        lo, hi = min(vals), max(vals)
        reciprocal = I(lo.numerator//lo.denominator,
                       -((-hi.numerator)//hi.denominator), True)
        return self*reciprocal

    def __rtruediv__(self, other):
        return I(other)/self

    def __pow__(self, k):
        assert isinstance(k, int) and k >= 0
        out, base = I(1), self
        while k:
            if k & 1:
                out = out*base
            k //= 2
            if k:
                base = base*base
        return out

    def __abs__(self):
        return I(max(0, self.lo, -self.hi),
                 max(abs(self.lo), abs(self.hi)), True)


def logunit(x):
    x = Q(x)
    assert 1 <= x <= 2
    z = (x-1)/(x+1)
    count = 55
    partial = sum((2*z**(2*j+1)/(2*j+1) for j in range(count)), Q(0))
    tail = 2*z**(2*count+1)/((2*count+1)*(1-z*z))
    return I.bounds(partial, partial+tail)


LOG2 = logunit(Q(2))


@lru_cache(None)
def logrec(r):
    x, m = 1/Q(r), 0
    while x > 2:
        x /= 2
        m += 1
    while x < 1:
        x *= 2
        m -= 1
    return m*LOG2+logunit(x)


class D:
    """An interval value and an interval enclosure of its first derivative."""
    def __init__(self, val=0, deriv=0):
        if isinstance(val, D):
            self.val, self.deriv = val.val, val.deriv
        else:
            self.val, self.deriv = I(val), I(deriv)

    def __add__(self, other):
        other = D(other)
        return D(self.val+other.val, self.deriv+other.deriv)

    __radd__ = __add__

    def __neg__(self):
        return D(-self.val, -self.deriv)

    def __sub__(self, other):
        return self+-D(other)

    def __rsub__(self, other):
        return D(other)+-self

    def __mul__(self, other):
        other = D(other)
        return D(self.val*other.val,
                 self.deriv*other.val+self.val*other.deriv)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = D(other)
        return D(self.val/other.val,
                 (self.deriv*other.val-self.val*other.deriv)/other.val**2)

    def __rtruediv__(self, other):
        return D(other)/self

    def __pow__(self, k):
        assert isinstance(k, int) and k >= 0
        return D(self.val**k, k*self.val**(k-1)*self.deriv) if k else D(1)


def neglog(x):
    x = D(x)
    assert x.val.lo > 0
    low = logrec(Q(x.val.hi, S))
    high = logrec(Q(x.val.lo, S))
    return D(I(low.lo, high.hi, True), -x.deriv/x.val)


def positive_series(x, offset):
    """exp(x) if offset=0; (exp(x)-1)/x if offset=1."""
    x = D(x)
    assert x.val.lo >= 0 and x.val.hi < 8*S
    value = I(0)
    deriv = I(0)
    power = I(1)
    previous = I(1)
    for j in range(ORDER+1):
        value += power/factorial(j+offset)
        if j:
            deriv += j*previous/factorial(j+offset)
        previous, power = power, power*x.val
    M = Q(x.val.hi, S)
    # All omitted terms are nonnegative; e^M < e^8 < 3^8.
    tail = I(3**8*M**(ORDER+1)/factorial(ORDER+1+offset))
    derivative_tail = I(3**8*M**ORDER/factorial(ORDER+offset))
    value = I(value.lo, value.hi+tail.hi, True)
    deriv = I(deriv.lo, deriv.hi+derivative_tail.hi, True)
    return D(value, deriv*x.deriv)


def low_test(q):
    q = D(q)
    denominator = 4+LINEAR*q+QUADRATIC*q*q
    r = q*q/denominator
    x = q*neglog(r)
    ex = positive_series(x, 0)
    c2 = 1+q*(1-q)/2
    b2 = 2+q*(1-q)**2/2
    return c2*ex+b2+1/ex-denominator


def high_test(q):
    q = D(q)
    r = q*q/(4+LINEAR*q+QUADRATIC*q*q)
    alpha = 1-q
    A = neglog(1-r)
    B = neglog(r)
    terms = [
        ((1-r)*(1+r+2*r*r), -(1-r)*(r+r*r), A),
        ((1-r)*r*(2+r+r*r), -(1-r)*r*(1+r), A+B),
        (-2*r*(1-r*r), r*(1-r*r), B),
        (r**4, D(0), 2*B),
    ]
    out = D(0)
    for a, b, lam in terms:
        x = alpha*lam
        out += a*lam*positive_series(x, 1)+b*positive_series(x, 0)
    # This is exactly (r+r^(1+q))^2 H'_(1+q)(r)/(1-q),
    # continuously extended to q=1, with no division by 1-q.
    return out


def certify(test, lo, hi, name):
    stack = [(Q(lo), Q(hi), 0)]
    accepted = 0
    margin = None
    max_depth = 0
    while stack:
        a, b, depth = stack.pop()
        mid, radius = (a+b)/2, (b-a)/2
        try:
            interval = test(D(I.bounds(a, b), 1))
            centre = test(D(mid, 0))
            upper = centre.val+abs(interval.deriv)*radius
        except AssertionError:
            upper = I(1)
        if upper.hi < 0:
            accepted += 1
            max_depth = max(max_depth, depth)
            margin = -upper.hi if margin is None else min(margin, -upper.hi)
        else:
            assert depth < 22, (name, a, b, Q(upper.hi, S))
            stack.append((mid, b, depth+1))
            stack.append((a, mid, depth+1))
    print(name, 'PASS', 'cells', accepted, 'depth', max_depth,
          'strict margin', Q(margin, S), flush=True)


if __name__ == '__main__':
    certify(low_test, Q(1,20), Q(1,2), 'low-q K upper bound')
    certify(high_test, Q(1,2), Q(1), 'high-q normalized derivative')
    print('PASS: r_p > 2q^2/(8+9q+23q^2) for every q in [1/20,1).')

