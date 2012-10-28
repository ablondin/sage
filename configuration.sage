# -*- coding: utf-8 -*-
#*****************************************************************************
#       Copyright (C) 2011 Alexandre Blondin Masse <ablondin@uqac.ca>,
#                          Ariane Garon <ariane.garon@gmail.com>,
#                          Sebastien Labbe <slabqc@gmail.com>
#
#  Distributed under the terms of the GNU General Public License
#  version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************
r"""
Class Configuration used to study the enumeration
of double square tiles.

AUTHORS:
- Alexandre Blondin Masse
- Sebastien Labbe
"""
class Configuration(object):
    r"""
    Configuration object.

    INPUT:

    - ``A`` - a word path having a 4-letter alphabet parent 
    - ``B`` - a word path having a 4-letter alphabet parent 
    - ``d1`` - integer
    - ``d2`` - integer

    or

    - ``A`` - a word path having a 4-letter alphabet parent 
    - ``B`` - a word path having a 4-letter alphabet parent 
    - ``d1`` - integer
    - ``d2`` - integer
    - ``bar`` - the bar operator

    or

    - ``ds`` - a double square

    or

    - ``l0`` - the length of w0
    - ``l1`` - the length of w1
    - ``l2`` - the length of w2
    - ``l3`` - the length of w3

    or

    - ``t`` - a tuple of 4 elements (w0,w1,w2,w3) or 8 elements
      (w0,w1,w2,w3,w4,w5,w6,w7)
    - ``bar`` - the bar operator

    EXAMPLES:

    From a configuration (A, B, d1, d2)::

        sage: P = WordPaths('abAB')
        sage: Configuration(P('aba'), P('bAb'),1,1)
        Configuration(w0, w1, w2, w3):            
          w0 = Path: a                              w4 = Path: A
          w1 = Path: ba                             w5 = Path: BA
          w2 = Path: b                              w6 = Path: B
          w3 = Path: Ab                             w7 = Path: aB
        (|w0|, |w1|, |w2|, |w3|) = (1, 2, 1, 2)   
        (n0, n1, n2, n3)         = (0, 1, 0, 1)  

    From a double square::

        sage: fibo = words.fibonacci_tile
        sage: Configuration(fibo(1))
        Configuration(w0, w1, w2, w3):            
          w0 = Path: 3                              w4 = Path: 1
          w1 = Path: 03                             w5 = Path: 21
          w2 = Path: 0                              w6 = Path: 2
          w3 = Path: 10                             w7 = Path: 32
        (|w0|, |w1|, |w2|, |w3|) = (1, 2, 1, 2)   
        (n0, n1, n2, n3)         = (0, 1, 0, 1)   
        sage: Configuration(fibo(2))
        Configuration(w0, w1, w2, w3):            
          w0 = Path: 30323                          w4 = Path: 12101
          w1 = Path: 21232303                       w5 = Path: 03010121
          w2 = Path: 23212                          w6 = Path: 01030
          w3 = Path: 10121232                       w7 = Path: 32303010
        (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
        (n0, n1, n2, n3)         = (0, 0, 0, 0)   

      From four integers::

        sage: c = Configuration(4,2,4,2)
        sage: c
        Configuration(w0, w1, w2, w3):            
          w0 = Path: DCab                           w4 = Path: cdBA
          w1 = Path: DC                             w5 = Path: cd
          w2 = Path: BADC                           w6 = Path: abcd
          w3 = Path: BA                             w7 = Path: ab
        (|w0|, |w1|, |w2|, |w3|) = (4, 2, 4, 2)   
        (n0, n1, n2, n3)         = (1, 0, 1, 0)   
        sage: c = Configuration(2,1,2,1)
        sage: c
        Configuration(w0, w1, w2, w3):            
          w0 = Path: Ab                             w4 = Path: aB
          w1 = Path: A                              w5 = Path: a
          w2 = Path: BA                             w6 = Path: ba
          w3 = Path: B                              w7 = Path: b
        (|w0|, |w1|, |w2|, |w3|) = (2, 1, 2, 1)   
        (n0, n1, n2, n3)         = (1, 0, 1, 0)   

    """
    def __init__(self, *args):
        r"""
        Constructor.

        See :Configuration: for documentation.

        EXAMPLES::

            sage: Configuration(2,2,2,2)
            Configuration(w0, w1, w2, w3):            
              w0 = Path: ab                             w4 = Path: ab
              w1 = Path: BA                             w5 = Path: BA
              w2 = Path: ab                             w6 = Path: ab
              w3 = Path: BA                             w7 = Path: BA
            (|w0|, |w1|, |w2|, |w3|) = (2, 2, 2, 2)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
        """
        if len(args) == 4 and all(isinstance(a, (int, Integer)) for a in args):
            msg = "les delais doivent etre strictement positifs"
            assert args[0] + args[2] > 0 and args[1] + args[3] > 0, msg 
            ((w0,w1,w2,w3,w4,w5,w6,w7), self.bar) = overlap_moi_ca(*args)
            self._w = (w0,w1,w2,w3,w4,w5,w6,w7)
            alphabet = self.A.parent().alphabet()
            if not hasattr(alphabet,'cardinality')\
            or alphabet.cardinality() != 4:
                pass
        elif len(args) == 4:
            A, B, d1, d2 = args
            w0 = A[:d1]
            w1 = A[d1:]
            w2 = B[:d2]
            w3 = B[d2:]
            self._w = (w0,w1,w2,w3)
            self.bar = None
        elif len(args) == 5:
            A, B, d1, d2, self.bar = args
            w0 = A[:d1]
            w1 = A[d1:]
            w2 = B[:d2]
            w3 = B[d2:]
            self._w = (w0,w1,w2,w3)
        elif len(args) == 2 :
            (self._w, self.bar) = args
            assert len(self._w) in (4, 8)
        elif len(args) == 1:
            ds = args[0]
            f = find_good_ds_factorisation(ds, delay='minimize')
            startA,endA,startX,endX = f
            demiper = ds.length()//2
            twice = ds * ds
            w0 = twice[startA:startX]
            w1 = twice[startX:endA]
            w2 = twice[endA:endX]
            w3 = twice[endX:startA+demiper]
            w4 = twice[startA+demiper:startX+demiper]
            w5 = twice[startX+demiper:endA+demiper]
            w6 = twice[endA+demiper:endX+demiper]
            w7 = twice[endX+demiper:startA+demiper+demiper]
            self._w = (w0,w1,w2,w3,w4,w5,w6,w7)
            self.bar = None
        else:
            raise TypeError, "Configuration takes one argument (a double square)\
                  or four arguments (A, B, d1, d2) not %s."%len(args)

        if self.bar is None:
            if self.A.parent() != self.B.parent():
                raise ValueError, "A and B must have the same parent"
            alphabet = self.A.parent().alphabet()
            if not hasattr(alphabet,'cardinality')\
            or alphabet.cardinality() != 4:
                raise ValueError, "The parent of A must have a 4-letter\
                                   alphabet."
            e,n,w,s = alphabet
            self.bar = WordMorphism({e:w,w:e,n:s,s:n},codomain=self.A.parent())

        self.verify_definition()
        self.verify_conjecture()

    @lazy_attribute
    def hat(self):
        return lambda x:self.bar(x).reversal()

    @lazy_attribute
    def A(self):
        return self._w[0] * self._w[1]

    @lazy_attribute
    def B(self):
        return self._w[2] * self._w[3]

    @lazy_attribute
    def d1(self):
        return len(self._w[0])

    @lazy_attribute
    def d2(self):
        return len(self._w[2])

    @cached_method
    def __getitem__(self, i):
        r"""
        Return the factor w_i 

        This corresponds to the new definition of configuration (solution).

        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(1))
            sage: [c[i] for i in range(8)]
            [Path: 3, Path: 03, Path: 0, Path: 10, Path: 1, Path: 21, Path: 2,
             Path: 32]

        ::

            sage: c = Configuration(fibo(2))
            sage: [c[i] for i in range(8)]
            [Path: 30323, Path: 21232303, Path: 23212, Path: 10121232,
             Path: 12101, Path: 03010121, Path: 01030, Path: 32303010]
        """
        if 0 <= i < len(self._w):
            return self._w[i]
        elif i == 4:
            return self.hat(self.A)[:self.d1]
        elif i == 5:
            return self.hat(self.A)[self.d1:]
        elif i == 6:
            return self.hat(self.B)[:self.d2]
        elif i == 7:
            return self.hat(self.B)[self.d2:]
        else:
            raise ValueError, 'i (=%s) must be between 0 and 7.'%i

    w = __getitem__

    def __eq__(self, other):
        r"""
        Returns True if A, B, d1 and d2 are the same.

        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c == c
            True
            sage: c == c.conjugate()
            False
            sage: c == c.old_extend(2,4).old_shrink(2,4)
            True
        """
        return isinstance(other, Configuration) and\
               self.A == other.A and\
               self.B == other.B and\
               self.d1 == other.d1 and\
               self.d2 == other.d2 

    def __cmp__(self, other):
        self_bw = self.boundary_word()
        other_bw = other.boundary_word()
        if len(self_bw) != len(other_bw):
            return len(self_bw) - len(other_bw)
        else:
            return self_bw.__cmp__(other_bw)

    def __len__(self):
        return len(self.boundary_word())

    def __hash__(self):
        r"""
        hash

        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(2))
            sage: hash(c)
            368453717
            sage: hash(c.exchange(0).exchange(0))
            368453717
        """
        return hash((self.A,self.B,self.d1,self.d2))

    def verify_definition(self):
        r"""
        Checks that the solution verify the definition.
        """
        for i in range(4):
            msg = "wiwi+1 = hat(wi+4,wi+5) is broken for i=%s"%i
            assert self[i] * self[i+1] ==\
                   self.hat(self[i+4] * self[(i+5)%8]), msg

    def verify_conjecture(self):
        r"""
        Verification de conjectures.

        Si une conjecture nest pas verifie, une AssertionError est lancee.
        """
        #Verification que chapeau de uivi est facteur de w pour tout i
        bw = self.boundary_word()
        for i in range(8):
            factor = self.hat(self.u(i)*self.v(i))
            a = factor.is_factor(bw**2)
            assert a, "hat(uivi)=(%s) nest pas facteur du contour pour i=%s\
                       et pour %s"%(factor,i,self)
        #Verification que chapeau de viui est facteur de w pour tout i
        for i in range(8):
            factor = self.hat(self.v(i)*self.u(i))
            a = factor.is_factor(bw**2)
            assert a, "hat(viui)=(%s) nest pas facteur du contour pour i=%s\
                       et pour %s"%(factor,i,self)

    def alphabet(self):
        r"""
        Returns the python set of the letters that occurs in the boundary
        word.
        """
        return set(self.boundary_word())

    def boundary_word(self):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.boundary_word()
            Path: 3032321232303232121012123212101030101210...
        """
        return self.A*self.B*self.hat(self.A)*self.hat(self.B)

    def turning_number(self):
        r"""
        """
        boundary = self.boundary_word()
        boundary = boundary + boundary[:1]
        boundary = boundary.to_integer_word()
        turns = boundary.finite_differences(mod=4)
        ev = turns.evaluation_dict()
        return QQ(((ev[1] - ev[3]), 4))

    def __repr__(self):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
        """
        s = []
        s.append("Configuration(w0, w1, w2, w3):")
        s += ["  w%s = %s"%(i,repr(self[i])) for i in range(4)]
        t = [''] + ["  w%s = %s"%(i,repr(self[i])) for i in range(4,8)]
        s.append("(|w0|, |w1|, |w2|, |w3|) = %s"%(tuple(map(len,\
                 (self[i] for i in range(4)))),))
        s.append("(n0, n1, n2, n3)         = %s"%(tuple(self.n(i)\
                                                  for i in range(4)),))
        return tableau_de_col(s,t)

    def _latex_(self):
        r"""
        Returns a 8-tuple representing the configuration in Latex

        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: latex(c)
            (30323,21232303,23212,10121232,12101,03010121,01030,32303010)

        """
        from sage.misc.latex import LatexExpr
        w = [self[i].string_rep() for i in range(8)]
        w = map(lambda s: '\\' + '\\'.join(s), w)
        return LatexExpr("(%s)"%", ".join(w))

    @cached_method
    def u(self, i):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c.u(0)
            Path: 30323
            sage: c.v(0)
            Path: 21232303010
            sage: c[0]
            Path: 30323
        """
        p = self.hat(self[(i-3)%8]) * self[(i-1)%8]
        return p[:len(self[i])%len(p)]

    @cached_method
    def v(self, i):
        r"""
        EXAMPLES::
        """
        p = self.hat(self[(i-3)%8]) * self[(i-1)%8]
        return p[len(self[i])%len(p):]

    @cached_method
    def n(self, i):
        r"""
        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: [c.n(i) for i in range(8)] 
            [0, 0, 0, 0, 0, 0, 0, 0]

        """
        p = self.hat(self[(i-3)%8]) * self[(i-1)%8]
        return len(self[i]) // len(p)

    @cached_method
    def d(self, i):
        r"""
        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0) 
            sage: [c.d(i) for i in range(8)]
            [16, 10, 16, 10, 16, 10, 16, 10]
        """
        return len(self.w((i - 1) % 8)) + len(self.w((i + 1) % 8))

    def thickness(self):
        r"""
        Returns the parameters (i,j) of the configuration.
        Note that thickness depends directly on the configuration,
        not only on the associated double square
        (see Words2009 article)

        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(2))
            sage: c.thickness()
            (0, 0)
            sage: c.old_extend(2,3).thickness()
            (2, 3)

            sage: c = Configuration(words.fibonacci_tile(2)).old_extend(2,3)
            sage: for i in range(8):
            ...     c = c.conjugate()
            ...     print c.thickness()
            (0, 0)
            (3, 2)
            (0, 0)
            (2, 3)
            (0, 0)
            (3, 2)
            (0, 0)
            (2, 3)
        """
        i = int((len(self.A) - self.d1) / (self.d1 + self.d2))
        j = int((len(self.B) - self.d2) / (self.d1 + self.d2))
        return (i,j)

    def permute(self):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.permute()
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 23212                          w4 = Path: 01030
              w1 = Path: 10121232                       w5 = Path: 32303010
              w2 = Path: 12101                          w6 = Path: 30323
              w3 = Path: 03010121                       w7 = Path: 21232303
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.permute().permute().permute().permute() == c
            True
        """
        return Configuration(self.B, self.hat(self.A),\
                             self.d2, self.d1, self.bar)

    def conjugate(self):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.conjugate()
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 21232303                       w4 = Path: 03010121
              w1 = Path: 23212                          w5 = Path: 01030
              w2 = Path: 10121232                       w6 = Path: 32303010
              w3 = Path: 12101                          w7 = Path: 30323
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.conjugate().conjugate().conjugate().conjugate().
                  conjugate().conjugate().conjugate().conjugate() == c
            True
            sage: c.conjugate().conjugate().conjugate().conjugate() == c
            False
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = tuple(self.w(i) for i in range(8))
        return Configuration((w1,w2,w3,w4,w5,w6,w7,w0), self.bar)

    def conjugates(self):
        conjugates = [self]
        for i in range(7):
            conjugates.append(conjugates[-1].conjugate())
        return conjugates

    def rotate(self):
        freeman = self.boundary_word().parent()
        a,b,A,B = freeman.alphabet()
        rho = WordMorphism({a:b, b:A, A:B, B:a}, codomain=freeman)
        return Configuration(tuple(rho(self.w(i)) for i in range(8)), self.bar)

    def reflect(self):
        freeman = self.boundary_word().parent()
        a,b,A,B = freeman.alphabet()
        sigma = WordMorphism({a:a, b:B, A:A, B:b}, codomain=freeman)
        return Configuration(tuple(sigma(self.w(i)) for i in range(8)), self.bar)

    def isometric_representant(self):
        if self.turning_number() != 1:
            candidates = [self.reverse()]
        else:
            candidates = [self]
        for i in range(3):
            candidates.append(candidates[-1].rotate())
        if self.turning_number() != 1:
            candidates.append(self.reflect())
        else:
            candidates.append(self.reverse().reflect())
        for i in range(3):
            candidates.append(candidates[-1].rotate())
        return min(map(lambda c: min(c.conjugates()), candidates))
            

    def reverse(self):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.reverse()
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 23212101                       w4 = Path: 01030323
              w1 = Path: 21232                          w5 = Path: 03010
              w2 = Path: 30323212                       w6 = Path: 12101030
              w3 = Path: 32303                          w7 = Path: 10121
            (|w0|, |w1|, |w2|, |w3|) = (8, 5, 8, 5)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.reverse().reverse() == c
            True
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = map(self.w, range(8))
        return Configuration((self.hat(w7),self.hat(w6),self.hat(w5),\
                              self.hat(w4),self.hat(w3),self.hat(w2),\
                              self.hat(w1),self.hat(w0)), self.bar)

    def extend(self, i):
        r"""
        EXAMPLES::
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = [self[j] for j in range(i % 8, 8) +\
                                     range(0, i % 8)]
        w = (w0*w1*self.hat(w3),w1,w2,w3,w4*w5*self.hat(w7),w5,w6,w7)
        w = tuple(w[j] for j in range(-i % 8, 8) + range(0, -i % 8))
        return Configuration(w, self.bar)

    def old_extend(self, k, l):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.old_extend(2,2)
            Configuration(w0, w1, w2, w3):                       
              w0 = Path: 3032321232303010303232123230301030323
              w1 = Path: 21232303
              w2 = Path: 2321210121232303232121012123230323212
              w3 = Path: 10121232
              w4 = Path: 1210103010121232121010301012123212101
              w5 = Path: 03010121
              w6 = Path: 0103032303010121010303230301012101030
              w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (37, 8, 37, 8)            
            (n0, n1, n2, n3)         = (2, 0, 2, 0)              
            sage: c.old_extend(3,5).old_shrink(3,5) == c
            True
        """
        c = self
        for _ in range(k): c = c.extend(0)
        for _ in range(l): c = c.extend(2)
        return c

    def shrink(self, i):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: d = c.old_extend(2,5)
            sage: d
            Configuration(w0, w1, w2, w3):                             
              w0 = Path: 3032321232303010303232123230301030323
              w1 = Path: 21232303
              w2 = Path: 2321210121232303232121012123230323212101...
              w3 = Path: 10121232
              w4 = Path: 1210103010121232121010301012123212101
              w5 = Path: 03010121
              w6 = Path: 0103032303010121010303230301012101030323...
              w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (37, 8, 85, 8)                  
            (n0, n1, n2, n3)         = (2, 0, 5, 0)                    
            sage: c == d.old_shrink(2,5) == d.old_shrink(1,3).old_shrink(1,2)
            True
        """
        d = len(self[(i-1)%8]) + len(self[(i+1)%8])
        if len(self[i]) > d:
            (w0,w1,w2,w3,w4,w5,w6,w7) = [self[j] for j in range(i % 8, 8)\
                                         + range(0, i % 8)]
            w = (w0[d:],w1,w2,w3,w4[d:],w5,w6,w7)
            w = tuple(w[j] for j in range(-i % 8, 8) + range(0, -i % 8))
            return Configuration(w, self.bar)
        else:
            raise ValueError, 'shrink cannot be applied on index %s\
                               of the following configuration\n%s'%(i,self)
    
    def old_shrink(self, k, l):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: d = c.old_extend(2,5)
            sage: d
            Configuration(w0, w1, w2, w3):                             
              w0 = Path: 3032321232303010303232123230301030323           
              w1 = Path: 21232303                                        
              w2 = Path: 2321210121232303232121012123230323212101...     
              w3 = Path: 10121232                                        
              w4 = Path: 1210103010121232121010301012123212101
              w5 = Path: 03010121
              w6 = Path: 0103032303010121010303230301012101030323...
              w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (37, 8, 85, 8)                  
            (n0, n1, n2, n3)         = (2, 0, 5, 0)                    
            sage: c == d.old_shrink(2,5) == d.old_shrink(1,3).old_shrink(1,2)
            True
        """
        c = self
        for _ in range(k): c = c.shrink(0)
        for _ in range(l): c = c.shrink(2)
        return c
    
    def exchange(self, i):
        r"""
        EXAMPLES::
        """
        (w0,w1,w2,w3,w4,w5,w6,w7) = [self[j] for j in range(i % 8, 8)\
                                     + range(0, i % 8)]
        indexes = range(i % 8 + 1, 8, 2) + range((1 + i) % 2, i % 8, 2)
        (n1,n3,n5,n7) = [self.n(j) for j in indexes]
        (u1,u3,u5,u7) = [self.u(j) for j in indexes]
        (v1,v3,v5,v7) = [self.v(j) for j in indexes]
        w = (self.hat(w4),(v1*u1)**n1*v1,self.hat(w6),\
             (v3*u3)**n3*v3,self.hat(w0),(v5*u5)**n5*v5,\
             self.hat(w2),(v7*u7)**n7*v7)
        w = tuple(w[j] for j in range(-i % 8, 8) + range(0, -i % 8))
        return Configuration(w, self.bar)

    def old_exchange(self):
        r"""
        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(2))
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.old_exchange()
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 32303                          w4 = Path: 10121
              w1 = Path: 23                             w5 = Path: 01
              w2 = Path: 21232                          w6 = Path: 03010
              w3 = Path: 12                             w7 = Path: 30
            (|w0|, |w1|, |w2|, |w3|) = (5, 2, 5, 2)   
            (n0, n1, n2, n3)         = (1, 0, 1, 0)   
            sage: c.old_exchange().old_exchange()
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 30323                          w4 = Path: 12101
              w1 = Path: 21232303                       w5 = Path: 03010121
              w2 = Path: 23212                          w6 = Path: 01030
              w3 = Path: 10121232                       w7 = Path: 32303010
            (|w0|, |w1|, |w2|, |w3|) = (5, 8, 5, 8)   
            (n0, n1, n2, n3)         = (0, 0, 0, 0)   
            sage: c.old_exchange().old_exchange() == c
            True
        """
        #shorter names
        A = self.A
        B = self.B
        d1 = self.d1
        d2 = self.d2
        d = d1 + d2

        A1 = A[:d1]
        A2 = A[-d1:]
        B1 = B[:d2]
        B2 = B[-d2:]
        r = (B2+A1)[:(len(A)-d1)%d]
        s = (B1+A2)[:d-len(r)]
        m = (self.hat(A1)+B1)[:(len(B)-d2)%d]
        n = (self.hat(A2)+B2)[:d-len(m)]
        i = int((len(A) - d1 - len(r))) / d
        j = int((len(B) - d2 - len(m))) / d

        Ap = A2 + (s + r) ** i + s
        Bp = B2 + (n + m) ** j + n

        return Configuration(Ap, Bp, d1, d2, self.bar)

    def old_mini_shrink(self):
        r"""
        Reduces the configuration in the case ``(|A| - d_1) \bmod (d_1 + d_2) =
        0`` and ``(|B| - d_2) \bmod (d_1 + d_2) \neq 0`` by using the smaller
        periodicity in the overlaps.

        EXAMPLES::

            sage: FREEMAN = WordPaths('abAB')
            sage: c = Configuration(FREEMAN('abAba'), FREEMAN('bAbAb'), 1, 3); c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: a                              w4 = Path: A
              w1 = Path: bAba                           w5 = Path: BaBA
              w2 = Path: bAb                            w6 = Path: BaB
              w3 = Path: Ab                             w7 = Path: aB
            (|w0|, |w1|, |w2|, |w3|) = (1, 4, 3, 2)   
            (n0, n1, n2, n3)         = (0, 1, 0, 0)   
            sage: c.old_mini_shrink()
            Configuration(w0, w1, w2, w3):            
              w0 = Path: a                              w4 = Path: A
              w1 = Path: ba                             w5 = Path: BA
              w2 = Path: b                              w6 = Path: B
              w3 = Path: Ab                             w7 = Path: aB
            (|w0|, |w1|, |w2|, |w3|) = (1, 2, 1, 2)   
            (n0, n1, n2, n3)         = (0, 1, 0, 1)   
        """
        (A,B,d1,d2) = (self.A,self.B,self.d1,self.d2)
        d = d1 + d2
        r = A[d1:d1+(len(A)-d1)%d]
        m = B[d2:d2+(len(B)-d2)%d]

        if len(r) == 0 and len(m) != 0:
            k = gcd(len(m), d - len(m))
            A1 = A[:d1]
            B1 = B[:d2]
            
            
            w = A1[:len(A1)%k]
            z = B1[:len(B1)%k]
            Ap = w + z + w
            Bp = z + m
            return Configuration(Ap, Bp, d1, len(z), self.bar)
        else:
            return Configuration(A, B, d1, d2, self.bar)

    def mini_shrink(self, i):
        r"""
        Reduces the configuration in the case ``|u_1| = 0``
        and ``|u_3| \neq 0`` by using the smaller
        periodicity in the overlaps.

        EXAMPLES::

            sage: c = Configuration(FREEMAN('abAba'), FREEMAN('bAbAb'), 1, 3).
            conjugate()
            sage: c
            Configuration(w0, w1, w2, w3):            
              w0 = Path: bAba                           w4 = Path: BaBA
              w1 = Path: bAb                            w5 = Path: BaB
              w2 = Path: Ab                             w6 = Path: aB
              w3 = Path: A                              w7 = Path: a
            (|w0|, |w1|, |w2|, |w3|) = (4, 3, 2, 1)   
            (n0, n1, n2, n3)         = (1, 0, 0, 0) 
            sage: c.u(0)
            Path: 
            sage: c.mini_shrink(2)
            Traceback (most recent call last):
            ...
            NotImplementedError: On ne peut enlever g(=2) a w_2, car |w_2|<=g
        """
        assert self.u(0).is_empty(), "u0 doit etre vide"
        w = w0,w1,w2,w3 = map(self.w, range(4))
        g = gcd(len(w2), len(w1) + len(w3))
        if g < len(w[i]):
            w[i] = w[i][g:]
        else:
            msg = "On ne peut enlever g(=%s) a w_%s, car |w_%s|<=g"%(g,i,i)
            raise NotImplementedError, msg
        return Configuration(w, self.bar)

    def l_shrink(self, i):
        r"""
        TODO
        """
        #assert len(self.w(i)) == self.d(i), 'on doit avoir w_%s == d_%s'%(i,i)
        w = [self.w((j+i)%8) for j in range(8)]
        d = [self.d((j+i)%8) for j in range(8)]
        g = gcd(len(w[2]),d[2]) 
        wp = [w[0][g:], w[1][g:], w[2], w[3], w[4][g:], w[5][g:], w[6], w[7]]
        return Configuration(tuple(wp[(j-i)%8] for j in range(8)), self.bar)

    def r_shrink(self, i):
        r"""
        TODO
        """
        #assert len(self.w(i)) == self.d(i), 'on doit avoir w_%s == d_%s'%(i,i)
        w = [self.w((j+i)%8) for j in range(8)]
        d = [self.d((j+i)%8) for j in range(8)]
        g = gcd(len(w[2]),d[2]) 
        wp = [w[0][:-g], w[1], w[2], w[3][:-g], w[4][:-g], w[5], w[6], w[7][:-g]]
        print tuple(wp[(j-i)%8] for j in range(8))
        return Configuration(tuple(wp[(j-i)%8] for j in range(8)), self.bar)

    def l_extend(self, i):
        r"""
        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(1))
            sage: c.l_extend(0)
            Traceback (most recent call last):
            ...
            AssertionError: on doit avoir w_0 == d_0
            sage: c.l_extend(1)
            Configuration(w0, w1, w2, w3):            
              w0 = Path: 3                              w4 = Path: 1
              w1 = Path: 0103                           w5 = Path: 2321
              w2 = Path: 010                            w6 = Path: 232
              w3 = Path: 10                             w7 = Path: 32
              (|w0|, |w1|, |w2|, |w3|) = (1, 4, 3, 2)   
              (n0, n1, n2, n3)         = (0, 1, 0, 0) 
        """
        assert len(self.w(i)) == self.d(i), 'on doit avoir w_%s == d_%s'%(i,i)
        w = [self.w((j+i)%8) for j in range(8)]
        d = [self.d((j+i)%8) for j in range(8)]
        g = gcd(len(w[2]),d[2]) 
        p = (w[1] + w[2] + w[3])[:g]
        q = (w[5] + w[6] + w[7])[:g]
        wp = [p + w[0], p + w[1], w[2], w[3], q + w[4], q + w[5], w[6], w[7]]
        return Configuration(tuple(wp[(j-i)%8] for j in range(8)), self.bar)

    def r_extend(self, i):
        r"""
        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(1))
            sage: c.r_extend(0)
            Traceback (most recent call last):
            ...
            AssertionError: wiwi+1 = hat(wi+4,wi+5) is broken for i=0
            sage: c.r_extend(1)
            Configuration(w0, w1, w2, w3):            
                w0 = Path: 323                            w4 = Path: 101
                w1 = Path: 0323                           w5 = Path: 2101
                w2 = Path: 0                              w6 = Path: 2
                w3 = Path: 10                             w7 = Path: 32
                (|w0|, |w1|, |w2|, |w3|) = (3, 4, 1, 2)   
                (n0, n1, n2, n3)         = (0, 1, 0, 0) 
        """
        #assert len(self.w(i)) == self.d(i), 'on doit avoir w_%s == d_%s'%(i,i)
        w = [self.w((j+i)%8) for j in range(8)]
        d = [self.d((j+i)%8) for j in range(8)]
        g = gcd(len(w[2]),d[2]) 
        p = (w[1] + w[2] + w[3])[:g]
        q = (w[5] + w[6] + w[7])[:g]
        wp = [w[0] + q, w[1], w[2], w[3] + p, w[4] + p, w[5], w[6], w[7] + q]
        return Configuration(tuple(wp[(j-i)%8] for j in range(8)), self.bar)

    def reduction(self, iteration=1, verbose=True):
        r"""
        Reduces the current configuration if it is possible

        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(3))
            sage: c.reduction(7)[1]
            exchange(0) applied
            shrink(0) applied
            shrink(2) applied
            exchange(0) applied
            shrink(0) applied
            shrink(2) applied
            not reducible, this is a morphic pentamino
            ['\\EXCHANGE_0', '\\SHRINK_0', '\\SHRINK_2', '\\EXCHANGE_0',\
             '\\SHRINK_0', '\\SHRINK_2']
        """
        # We start with the recursive call
        if iteration >= 2:
            conf,op = self.reduction(1, verbose=verbose)
            conf2,op2 = conf.reduction(iteration - 1, verbose)
            return (conf2, op + op2)

        # Now we check if SHRINK may be applied
        for i in range(8):
            if len(self.w(i)) > self.d(i):
                if verbose:
                    print 'shrink(%s) applied'%i
                return (self.shrink(i), ['\\SHRINK_{%s}'%i])

        # Then we verify if we have a morphic pentamino
        if (len(self.u(0)) == 0 and len(self.u(2)) == 0) or\
           (len(self.u(1)) == 0 and len(self.u(3)) == 0):
            if verbose:
                print 'not reducible, this is a morphic pentamino'
            return (self, [])

        # Otherwise, we try with EXCHANGE
        for i in range(2):
            if len(self.v((i + 1) % 8)) + len(self.v((i + 3) % 8)) <\
               len(self.u((i + 1) % 8)) + len(self.u((i + 3) % 8)):
                if verbose:
                    print 'exchange(%s) applied'%i
                return (self.exchange(i), ['\\EXCHANGE_{%s}'%i])

        raise ValueError, 'case not treated yet !!!\n%s'%self

    def old_reduction(self, iteration=1, verbose=True):
        r"""
        Reduces the current configuration to a smaller one, if possible
        """
        if iteration >= 2:
            return self.old_reduction(1, verbose=verbose).old_reduction\
                   (iteration - 1, verbose)
        c = self.find_thickest_equivalent_configuration()
        (A,B,d1,d2) = (c.A,c.B,c.d1,c.d2)
        d = d1 + d2
        dp = len(A) + len(B) - d
        r = (len(A) - d1) % d
        m = (len(B) - d2) % d
        (i,j) = c.thickness()
        if r == 0: i = max(0, i - 1)
        if m == 0: j = max(0, j - 1)
        if i >= 1 or j >= 1:
            rep = c.shrink(i,j)
            if verbose: 
                print 'shrink applied with parameters', i, 'and', j
                rep = (rep, "$\\SHRINK(%s, %s)$"%(i,j))
        elif r == 0 and m == 0:
            rep = None
            if verbose: 
                print 'not reductible : r and m both empty'
                rep = (rep, None)
        elif d == (len(A) + len(B)) / 2:
            rep = c
            if verbose: 
                print 'not reductible : case d = d\''
                rep = (rep, None)
        elif r + m > d:
            rep = c.exchange()
            if verbose: 
                print 'exchange applied'
                rep = (rep, "\\EXCHANGE")
        elif d2 % dp + d1 % dp > dp:
            rep = c.conjugate().exchange()
            if verbose: 
                print 'exchange applied to conjugate'
                rep = (rep, "$\\EXCHANGE\\circ\\CONJUGATE$")
        elif r == 0:
            rep = c.old_mini_shrink()
            if verbose: 
                print 'mini shrink applied'
                rep = (rep, "$\\MSHRINK$")
        elif m == 0:
            rep = c.conjugate().conjugate().old_mini_shrink()
            if verbose: 
                print 'mini shrink applied to conjugate'
                rep = (rep, "$\\MSHRINK\\circ\\CONJUGATE\\circ\\CONJUGATE$")
        else:
            raise Exception, "missing case in reduction."
        return rep

    def find_thickest_equivalent_configuration(self):
        r"""
        Returns the thickest configuration equivalent to this
        one.
        More precisely, using the operators reverse() and conjugate(),
        it chooses the configuration having thickness (i,j) such that
        i + j is maximized

        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(1)).conjugate().\
            old_extend(1,2).conjugate()
            sage: c.thickness()
            (0, 0)
            sage: c.find_thickest_equivalent_configuration().thickness()
            (2, 1)
        """
        candidates = [self, self.conjugate(), self.reverse(),\
                      self.reverse().conjugate()]
        distinct_thicknesses = set(map(lambda c:c.thickness(), candidates))
        if len(distinct_thicknesses) > 2:
            print 'this seems to be a counter-example to the constant thickness\
                   hypothesis.'
            print distinct_thicknesses
        return max(candidates, key=lambda c:sum(c.thickness()))

    def factorization_points(self):
        r"""
        Returns the eight factorization points of this configuration
        """
        return [0, self.d1, len(self.A), len(self.A)+self.d2,\
                len(self.A)+len(self.B), len(self.A)+len(self.B)+self.d1,\
                2*len(self.A)+len(self.B), 2*len(self.A)+len(self.B)+self.d2]

    def plot(self, pathoptions=dict(rgbcolor='black',thickness=3), 
         fill=True, filloptions=dict(rgbcolor='black',alpha=0.2),
         startpoint=True, startoptions=dict(rgbcolor='black',pointsize=100), 
         endarrow=True, arrowoptions=dict(rgbcolor='black',arrowsize=5,width=3),
         gridlines=False, gridoptions=dict(),
         axes=False):
        r"""
        Returns a 2d Graphics illustrating the double square tile associated to
        this configuration including the factorizations points. The options are
        the same as for instances of WordPaths

        INPUT:

        - ``pathoptions`` - (dict,
          default:dict(rgbcolor='red',thickness=3)), options for the
          path drawing

        - ``fill`` - (boolean, default: True), if fill is True and if
          the path is closed, the inside is colored

        - ``filloptions`` - (dict,
          default:dict(rgbcolor='red',alpha=0.2)), ptions for the
          inside filling

        - ``startpoint`` - (boolean, default: True), draw the start point?
        
        - ``startoptions`` - (dict,
          default:dict(rgbcolor='red',pointsize=100)) options for the
          start point drawing

        - ``endarrow`` - (boolean, default: True), draw an arrow end at the end?

        - ``arrowoptions`` - (dict,
          default:dict(rgbcolor='red',arrowsize=20, width=3)) options
          for the end point arrow

        - ``gridlines``- (boolean, default: False), show gridlines?
        
        - ``gridoptions`` - (dict, default: {}), options for the gridlines
        
        - ``axes`` - (boolean, default: False), options for the axes


        EXAMPLES:

        The cross of area 5 together with its double square factorization
        points::

            sage: c = Configuration(words.fibonacci_tile(1))
            sage: c.plot()
        """
        path = self.boundary_word()
        points = list(path.points())
        points = [map(RR, x) for x in points]
        G = path.plot(pathoptions, fill, filloptions, startpoint, startoptions,\
                      endarrow, arrowoptions, gridlines, gridoptions)
        i = 0
        for p in self.factorization_points():
            if i % 2 == 0: G += point(points[p],\
                                      pointsize=startoptions['pointsize'],\
                                      rgbcolor="red")
            else: G += point(points[p],\
                             pointsize=startoptions['pointsize'],\
                             rgbcolor="blue")
            i += 1
        return G

    def tikz_trajectory(self, step=1):
        r"""
        Returns a tikz string describing the double square induced by
        this configuration together with its factorization points

        The factorization points respectively get the tikz attribute 'first'
        and 'second' so that when including it in a tikzpicture environment,
        it is possible to modify the way those points appear.

        EXAMPLES::

            sage: c = Configuration(words.fibonacci_tile(1))
            sage: c.tikz_trajectory()
            \filldraw[-triangle 45, very thick, draw=black, fill=black!10]
            (0.000, 0.000) -- (0.000, -1.00) -- (1.00, -1.00) -- (1.00, -2.00)
            -- (2.00, -2.00) -- (2.00, -1.00) -- (3.00, -1.00) -- (3.00, 0.000)
            -- (2.00, 0.000) -- (2.00, 1.00) -- (1.00, 1.00) -- (1.00, 0.000)
            -- (0.000, 0.000);
            \foreach \i in {(0.0000, 0.0000), (1.000, -2.000), (3.000, -1.000),
            (2.000, 1.000)}
              \node at \i[first] {};
            \foreach \i in {(0.0000, -1.000), (2.000, -2.000), (3.000, 0.0000),
            (1.000, 1.000)}
              \node at \i[second] {};
        """
        from sage.all import n
        from sage.misc.latex import LatexExpr
        f = lambda x: n(x,digits=3)
        step = n(step, digits=4)
        points = map(lambda (x,y):(x*step,y*step),\
                     list(self.boundary_word().points()))
        l = [str(tuple(map(f, pt))) for pt in points]
        s = '\\filldraw[-triangle 45, very thick, draw=black, fill=black!10] '\
            + ' -- '.join(l) + ';'
        [a1, b1, a2, b2, a3, b3, a4, b4] = self.factorization_points()
        #s += '\n\\foreach \\x / \\y in {'
        #for p in [a1, a2, a3]:
        #    x,y = points[p]
        #    s += '%s/%s, ' % (str(x), str(y))
        #x,y = points[a4]
        #s += '%s/%s' % (str(x), str(y))
        #s += '}'
        #s += '\n  \\node[first] at (\\x, \\y) {};\n'
        #s += '\n\\foreach \\x / \\y in {'
        #for p in [b1, b2, b3]:
        #    x,y = points[p]
        #    s += '%s/%s, ' % (str(x), str(y))
        #x,y = points[b4]
        #s += '%s/%s' % (str(x), str(y))
        #s += '}'
        #s += '\n  \\node[second] at (\\x, \\y) {};\n'
        return LatexExpr(s)

    def tikz_reduction(self, size=6, nbcolonnes=3):
        r"""
        INPUT:

        - ``nbcolonnes`` - le nombre de colonnes de l'affichage

        EXAMPLES::

            sage: fibo = words.fibonacci_tile
            sage: c = Configuration(fibo(1))
            sage: c.tikz_reduction()
            \node (q0) at (0, 0.000000000000000) {\begin{tikzpicture}
            \filldraw[-triangle 45, very thick, draw=black, fill=black!10]
            (0.000, 0.000) -- (0.000, -2.00) -- (2.00, -2.00) -- (2.00, -4.00)
            -- (4.00, -4.00) -- (4.00, -2.00) -- (6.00, -2.00) -- (6.00, 0.000)
            -- (4.00, 0.000) -- (4.00, 2.00) -- (2.00, 2.00) -- (2.00, 0.000)
            -- (0.000, 0.000);
            \foreach \i in {(0.0000, 0.0000), (2.000, -4.000), (6.000, -2.000),
            (4.000, 2.000)}
              \node at \i[first] {};
            \foreach \i in {(0.0000, -2.000), (4.000, -4.000), (6.000, 0.0000),
            (2.000, 2.000)}
              \node at \i[second] {};\end{tikzpicture}};
        """
        c = self
        s = ''
        i = 0
        fns = []
        while True:
            (w, h) = (c.width(), c.height())
            t = c.tikz_trajectory(step=size/max(w,h))
            x,y = serpent(i, nbcolonnes)
            s += '\\node (q%s) '%i
            s += 'at %s '% ( (x*1.5*size, y*1.5*size), )
            s += '{\\begin{tikzpicture}%s\\end{tikzpicture}};\n'%t
            c,func = c.reduction(iteration=1, verbose=False)
            i += 1
            if func == []: break
            fns.append(func[0])
        for j,func in zip(range(1, i), fns):
            rotate90 = "" if j%3==0 else ", rotate=90" 
            edge = "edge node[midway, rectangle, fill=white%s] "%rotate90
            edge += "{$%s$}"%func
            s += "\\path[->] (q%s) %s (q%s);\n"%(j-1, edge, j)
        from sage.misc.latex import LatexExpr
        return LatexExpr(s)

    def width(self):
        r"""
        Returns the width of this polyomino, i.e. the difference
        between its rightmost and leftmost coordinates
        """
        points = list(self.boundary_word().points())
        return max(map(lambda p:p[0], points)) - min(map(lambda p:p[0], points))

    def height(self):
        r"""
        Returns the width of this polyomino, i.e. the difference
        between its uppermost and lowermost coordinates
        """
        points = list(self.boundary_word().points())
        return max(map(lambda p:p[1], points)) - min(map(lambda p:p[1], points))

    def uv_non_simple_dict(self):
        r"""
        Retourne un dictionnaire de la forme 
        {'u': liste d'entiers, 'v': liste d'entiers}
        ou les liste d'entiers donne les indices i tels que ui ou vi se
        croisent.

        EXAMPLES::

            sage: c = Configuration(1,4,3,2)
            sage: c.uv_non_simple_dict()
            {'u': [], 'v': [0, 4]}
        """
        d = {}
        d['u'] = [i for i in range(8) if not self.u(i).is_simple() ]
        d['v'] = [i for i in range(8) if not self.v(i).is_simple() ]
        return d

    def uv_closed_dict(self):
        r"""
        Retourne un dictionnaire de la forme 
        {'u': liste d'entiers, 'v': liste d'entiers}
        ou les liste d'entiers donne les indices i tels que ui ou vi sont
        fermes (mais non vide!).

        EXAMPLES::

            sage: c = Configuration(1,4,3,2)
            sage: c.uv_closed_dict()
            {'u': [], 'v': []}

        ::

            sage: c = Configuration(words.fibonacci_tile(1))
            sage: d = c.conjugate().exchange(0)
            sage: d.uv_closed_dict()
            {'u': [], 'v': [0, 2, 4, 6]}
        """
        d = {}
        d['u'] = [i for i in range(8) if self.u(i).is_closed() and\
                  not self.u(i).is_empty() ]
        d['v'] = [i for i in range(8) if self.v(i).is_closed() and\
                  not self.v(i).is_empty() ]
        return d

    def latex_table(self):
        r"""
        Returns a Latex expression of a table containing
        the parameters A, B, X, Y, d1, d2, d, d'1, d'2,
        d', r, m, i and j of this configuration
        """
        from sage.misc.latex import LatexExpr
        remove_coma = lambda s:s.translate(None, ',')
        if_empty = lambda s:'\\varepsilon' if len(s) == 0 else s
        u = [if_empty(remove_coma(self.u(i).string_rep())) for i in range(4)]
        v = [if_empty(remove_coma(self.v(i).string_rep())) for i in range(4)]
        #uv = [self.u(i) for i in range(4)] + [self.v(i) for i in range(4)]
        x = self.uv_non_simple_dict()
        y = self.uv_closed_dict()
        secroisent = ', '.join(['u_%s'%i for i in x['u']] + ['v_%s'%i\
                     for i in x['v']])
        fermes = ', '.join(['u_%s'%i for i in y['u']] + ['v_%s'%i\
                     for i in y['v']])

        s = '\\begin{tabular}{|c|}\n\\hline\n\\\\\n'
        s += '\\begin{tikzpicture}\n'
        s += '  [first/.style={circle,draw=black,fill=gray, inner sep=0pt,\
             minimum size=3pt},\n'
        s += '   second/.style={rectangle,draw=black,fill=white, inner sep=0pt,\
             minimum size=3pt}]\n'
        s += self.tikz_trajectory(step=5.0/max(self.width(), self.height()))
        s += '\n\\end{tikzpicture} \\\\[1ex] \n\\hline\\\\\n'
        #for i in range(4):
        #    s += '$w_{%s} = %s$\\\\\n'%(i, remove_coma(self[i].string_rep()))
        s += '$(w_0,w_1,w_2,w_3) = (%s,%s,%s,%s)$ \\\\\n'%\
             tuple(len(self[i]) for i in range(4))
        s += '$u_0 = %s$\quad $u_1 = %s$\\\\$u_2 = %s$\quad $u_3 = %s$\\\\\n'%\
             tuple(u[i] for i in range(4))
        s += '$v_0 = %s$\quad $v_1 = %s$\\\\$v_2 = %s$\quad $v_3 = %s$\\\\\n'%\
             tuple(v[i] for i in range(4))
        s += '$(n_0,n_1,n_2,n_3) = (%s,%s,%s,%s)$ \\\\\n'%\
             tuple(self.n(i) for i in range(4))
        s += 'Turning number = %s\\\\\n'%self.turning_number()
        s += 'Self-avoiding = %s\\\\\n'%self.boundary_word().is_simple()
        s += 'Se croisent $=%s$\\\\\n'%secroisent
        s += 'Sont fermes $=%s$\\\\\n'%fermes
        s += '\\hline\n\\end{tabular}\n'
        return LatexExpr(s)

    ####################################
    # Fonctions sur les configurations #
    ####################################

    def is_factor_of_ui_or_vi(self, w):
        r"""
        Returns True if w is a factor of one of the ui or vi of self.
        """
        return any(w.is_factor(self.u(j)) or w.is_factor(self.v(j))\
                   for j in range(8))

    def uv_conservation_for_exchange_dict(self):
        r"""
        Retourne un dictionnaire qui donne les indices i tels que ui, vi et
        leurs chapeaux sont conserves par l'operateur exchange.

        OUTPUT:

            Un dictionaire de la forme 
            {'u': liste d'entiers, 'v': liste d'entiers, 'hatu': liste
            d'entiers, 'hatv': liste d'entiers}

        EXAMPLES::

            sage: c = Configuration(FREEMAN('abAba'), FREEMAN('bAbAb'),1,3).
            old_extend(1,1).conjugate().old_extend(1,2)
            sage: c.uv_conservation_for_exchange_dict()
            {'hatu': [0, 1, 2, 3, 4, 5, 6, 7], 'hatv': [0, 2, 4, 6],
            'u': [0, 1, 2, 3, 4, 5, 6, 7], 'v': [0, 1, 2, 3, 4, 5, 6, 7]}
        """
        e = self.exchange(0)
        hat = self.hat
        d = {}
        d['u'] = [i for i in range(8) if e.is_factor_of_ui_or_vi(self.u(i))]
        d['v'] = [i for i in range(8) if e.is_factor_of_ui_or_vi(self.v(i))]
        d['hatu'] = [i for i in range(8)\
                     if e.is_factor_of_ui_or_vi(hat(self.u(i)))]
        d['hatv'] = [i for i in range(8)\
                     if e.is_factor_of_ui_or_vi(hat(self.v(i)))]
        return d

    def checks_uv_conservation_for_exchange_operator(self, hat_allowed=True):
        r"""
        Checks which `u_i` and `v_i` of the given configuraiton
        are preserved when applying the exchange operator

        INPUT:

        - ``self`` - a configuration.
        - ``hat_allowed`` - a boolean (default: ``True``). If True, then cheks
          also if `\hat{u_i}` and `\hat{v_i}` occur.

        EXAMPLES::

            sage: c = Configuration(FREEMAN('abAba'), FREEMAN('bAbAb'),1,3).
            old_extend(1,1).conjugate().old_extend(1,2) 
            sage: c.checks_uv_conservation_for_exchange_operator() 
            ['u_0', 'u_1', 'u_2', 'u_3', 'u_4', 'u_5', 'u_6', 'u_7', 'v_0',\
             'v_1', 'v_2', 'v_3', 'v_4', 'v_5', 'v_6', 'v_7', '\\hat{u_0}',\
             '\\hat{u_1}', '\\hat{u_2}', '\\hat{u_3}', '\\hat{u_4}',\
             '\\hat{u_5}', '\\hat{u_6}', '\\hat{u_7}', '\\hat{v_0}',\
             '\\hat{v_2}', '\\hat{v_4}', '\\hat{v_6}']

        """
        uv = []
        d = self.uv_conservation_for_exchange_dict()
        uv.extend('u_%s'%j for j in d['u'])
        uv.extend('v_%s'%j for j in d['v'])
        if hat_allowed:
            uv.extend('\\hat{u_%s}'%j for j in d['hatu'])
            uv.extend('\\hat{v_%s}'%j for j in d['hatv'])
        return uv

###############################
# String manipulation helpers #
###############################

def tableau_de_col(col1, col2, espace=3):
    r"""
    EXAMPLES::

        sage: col1 = ['ab','asdfasdf', 'adf']
        sage: col2 = ['11', '1313', '131313', '1313']
        sage: print tableau_de_col(col1, col2)
        ab         11
        asdfasdf   1313
        adf        131313
                   1313
    """
    from itertools import izip_longest
    largeur = max(map(len, col1))
    it = izip_longest(col1,col2, fillvalue='')
    espace = ' '*espace
    L = [a.ljust(largeur) + espace + b for (a,b) in it]
    return '\n'.join(L)

##################################
# Finding a square factorization #
##################################

def find_square_factorisation(ds, factorisation=None, alternate=True):
    r"""
    Return a square factorisation of the double square ds, distinct from
    factorisation.

    INPUT:

    - ds - word, a tile
    - factorisation - tuple (optional), a known factorisation
    - alternate - bool (optional, default True), if True the search for the
      second factorisation is restricted to those who alternates with the
      first factorisation

    OUTPUT:

    tuple of four positions of a square factorisation

    EXAMPLES::

        sage: fibo = words.fibonacci_tile
        sage: find_square_factorisation(fibo(1))
        (0, 3, 6, 9)
        sage: find_square_factorisation(fibo(0))
        (0, 1, 2, 3)
        sage: find_square_factorisation(fibo(1))
        (0, 3, 6, 9)
        sage: find_square_factorisation(fibo(2))
        (0, 13, 26, 39)
        sage: find_square_factorisation(fibo(3))
        (0, 55, 110, 165)

    ::

        sage: f = find_square_factorisation(fibo(3));f
        (0, 55, 110, 165)
        sage: find_square_factorisation(fibo(3),f)
        (34, 89, 144, 199)
        sage: find_square_factorisation(fibo(3),f,False) #optional long
        (34, 89, 144, 199)

    ::

        sage: find_square_factorisation(christo_tile(4,5))
        (0, 7, 28, 35)
        sage: find_square_factorisation(christo_tile(4,5),_)
        (2, 27, 30, 55)

    ::

        sage: find_square_factorisation(Words('abcd')('aaaaaa'))
        Traceback (most recent call last):
        ...
        ValueError: pas de factorisation carree
        sage: find_square_factorisation(Words('abcd')('aaaaaa'),(1,2,3,4))
        Traceback (most recent call last):
        ...
        ValueError: pas de seconde factorisation carree

    """
    e,n,w,s = ds.parent().alphabet()
    bar = WordMorphism({e:w,w:e,n:s,s:n},codomain=ds.parent())
    hat = lambda x:bar(x).reversal()
     
    l = ds.length()
    demiper = l/2
    aucarre = ds * ds
     
    if factorisation and alternate:
        a,b,c,d = factorisation
        it = ((debutA,finA) for debutA in range(a+1,b)\
              for finA in range(b+1,a+demiper))
    else:
        it = ((debutA,finA) for debutA in range(demiper)\
              for finA in range(debutA+1,debutA+demiper+1) )
    
    for debutA,finA in it:
        new = (debutA,finA,(debutA+demiper)%l,(finA+demiper)%l)
        if factorisation and set(factorisation) == set(new):
            continue
        A = aucarre[debutA:finA]
        B = aucarre[finA:debutA+demiper]
        A2 = aucarre[debutA+demiper:finA+demiper]
        B2 = aucarre[finA+demiper:l+debutA]
        if A == hat(A2) and B == hat(B2):
            return new
    
    if factorisation is None:
        raise ValueError, 'pas de factorisation carree'
    else:
        raise ValueError, 'pas de seconde factorisation carree'

def find_good_ds_factorisation(ds, delay='minimize'):
    r"""
    Returns the factorizations such that d(A,X) + d(B,Y) is 
    minimized (or maximized) where d(x,y) is the distance between 
    the start of the block x and the start of the block y.

    INPUT:

    - ``ds`` - double square
    - ``delay`` - 'minimize' or 'maximize'

    OUTPUT:

    tuple of four integers (start of A, end of A, start of X, end of X)

    EXAMPLES::

        sage: fibo = words.fibonacci_tile
        sage: find_good_ds_factorisation(fibo(1))
        (2, 5, 3, 6)
        sage: find_good_ds_factorisation(fibo(2))
        (8, 21, 13, 26)
        sage: find_good_ds_factorisation(fibo(3))
        (34, 89, 55, 110)
        sage: find_good_ds_factorisation(fibo(1),delay='maximize')
        (0, 3, 2, 5)
        sage: find_good_ds_factorisation(fibo(2),delay='maximize')
        (0, 13, 8, 21)
        sage: find_good_ds_factorisation(fibo(3),delay='maximize') #optional long
        (0, 55, 34, 89)

    """
    f = find_square_factorisation(ds)
    g = find_square_factorisation(ds, f, alternate=True)
    #print f,g
    debutA,finA,debutAhat,finAhat = sorted(f)
    debutX,finX,debutXhat,finXhat = sorted(g)

    l = ds.length()
    demiper = l//2
    d1 = debutX - debutA
    d2 = finX - finA
    #print d1,d2,demiper

    if (delay == 'minimize' and d1 + d2 > demiper/2) or\
       (delay == 'maximize' and d1 + d2 < demiper/2):
        #change the factorisation: B becomes X, X becomes A
        debutX,finX,debutA,finA = finA,debutA+demiper,debutX,finX
        d1 = debutX - debutA
        d2 = finX - finA
    #print d1,d2,demiper
    return debutA,finA,debutX,finX

#-------------------------------#
# Enumeration of double squares #
#-------------------------------#

Freeman = words.fibonacci_tile(1).parent()
F_BAR   = WordMorphism({0:2, 1:3, 2:0, 3:1}, codomain=Freeman)
F_HAT   = lambda w: F_BAR(w.reversal())
F_RHO   = WordMorphism({0:1, 1:2, 2:3, 3:0}, codomain=Freeman)
F_SIGMA = WordMorphism({0:2, 1:1, 2:0, 3:3}, codomain=Freeman)

def isometric_paths(path):
    return [path, F_RHO(path), F_RHO(F_RHO(path)), F_RHO(F_RHO(F_RHO(path))),\
            F_SIGMA(path), F_RHO(F_SIGMA(path)), F_RHO(F_RHO(F_SIGMA(path))),\
            F_RHO(F_RHO(F_RHO(F_SIGMA(path))))]

def all_double_squares_iterator(max_iteration=None,\
                                max_length=None,\
                                verbose=False):
    if max_iteration is None:
        max_iteration = Infinity
    if max_length is None:
        max_length = Infinity
    import heapq
    queue = [Configuration(words.fibonacci_tile(1))]
    visited = set([])
    i = 0
    while queue and i < max_iteration:
        tile = heapq.heappop(queue).isometric_representant()
        if not tile in visited:
            if tile.boundary_word().is_simple():
                if verbose:
                    print i, tile.boundary_word()
                yield tile.boundary_word()
                i += 1
            visited |= set([tile])
            # Extend operator
            for j in range(4):
                t = tile.extend(j)
                if len(t.boundary_word()) <= max_length:
                    heapq.heappush(queue, t)
            # Extend operator
            for j in range(2):
                t = tile.exchange(j)
                if len(t) >= len(tile) and len(t.boundary_word()) <= max_length:
                    heapq.heappush(queue, t)
