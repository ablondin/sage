#*****************************************************************************
#       Copyright (C) 2011 Alexandre Blondin Masse <ablondin@uqac.ca>,
#
#  Distributed under the terms of the GNU General Public License
#  version 2 (GPLv2)
#
#  The full text of the GPLv2 is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

#----------#
# Constant #
#----------#

Freeman = WordPaths('abAB')
BAR     = lambda n: -n
F_BAR   = WordMorphism('a->A,b->B,A->a,B->b', codomain=Freeman)
F_HAT   = lambda w: F_BAR(w.reversal())
F_RHO   = WordMorphism('a->b,b->A,A->B,B->a', codomain=Freeman)
F_SIGMA = WordMorphism('a->A,b->b,A->a,B->B', codomain=Freeman)

def PHI(n):
    """
    This morphism is used to transform general words into
    words on the Freeman alphabet
    """
    if n == 1:
        return 'a'
    elif n == -1:
        return 'A'
    elif n == 2:
        return 'b'
    elif n == -2:
        return 'B'

#---------------------#
# Partition functions #
#---------------------#

def get_generic_partition(generic_expressions, delays):
    """
    Returns the equivalence classes induced by overlaping
    the given generic expressions.

    EXAMPLES::

        sage: get_generic_partition([[1,2,3,4,5],[5,4,3,2,1]], [0])
        [set([1,5]),set([2,4]),set([3])]
    """
    l = min(delays + [0])
    u = len(generic_expressions[0])
    i = 1
    while i < len(generic_expressions):
        u = max(u, len(generic_expressions[i]) + delays[i - 1])
        i += 1
    partition = []
    for i in range(l, u):
        if i in range(len(generic_expressions[0])):
            c = set([generic_expressions[0][i]])
        else:
            c = set([])
        for j in range(1, len(generic_expressions)):
            if i in range(delays[j - 1],\
                          len(generic_expressions[j]) + delays[j - 1]):
                c |= set([generic_expressions[j][i - delays[j - 1]]])
        i = 0
        partition = merge_group_in_partition(c, partition)
        partition = merge_group_in_partition(set(map(lambda n:-n, c)), partition)
    return partition

def merge_group_in_partition(group, partition):
    """
    Merges the given group to the partition
    If this group creates bigger classes, they are merged as well.

    EXAMPLES:
        sage: merge_group_in_partition([set([2,3]), [set([1,2]),set([3,4])])
        [set([1,2,3,4])]
        sage: merge_group_in_partition([set([1,2]), [set([1,2]),set([3,4])])
        [set([1,2]),set([3,4])]
        sage: merge_group_in_partition([set([5]), [set([1,2]),set([3,4])])
        [set([1,2]),set([3,4]),set([5])]
    """
    if len(group) == 0:
        return copy(partition)
    else:
        merged_partition = []
        merged_class = group.copy()
        for c in partition:
            if len(c & group) == 0:
                merged_partition.append(c)
            else:
                merged_class |= c
        merged_partition.append(merged_class)
        return merged_partition

def get_representants(partition):
    """
    Returns a dictionary giving, for each number its representant
    according to the given partition.
    For sake of simplicity, the integer close to 0 is chosen.

    EXAMPLES:
        sage: get_representants([set([1,2]),set([3,4,5])])
        {1 : 1, 2 : 1, 3 : 3, 4 : 3, 5 : 3}
    """    
    representants = {}
    for p in partition:
        r = min(p, key = lambda n:abs(n))
        for q in p:
            representants[q] = r
    return representants

def is_realizable(partition):
    """
    Returns True if no letter is equivalent to its complement

    EXAMPLES:
        sage: is_realizable([set([1,2]),set([3,4])])
        True
        sage: is_realizable([set([1,2,-1]),set([3,4])])
        False
    """

    realizable = True
    i = 0
    while realizable and i < len(partition):
        realizable = len(set(map(lambda x:abs(x), partition[i]))) ==\
                     len(partition[i])
        i += 1
    return realizable

#----------------------------#
# Solving equations on words #
#----------------------------#

def get_generic_words(lengths):
    """
    Returns the general solution of the given equation on words.

    EXAMPLES:
        sage: get_generic_words({'u' : 4})
        {'u' : [1,2,3,4]}
        sage: get_generic_words({'u' : 2, 'v' : 6})
        {'u' : [1,2], 'v' : [3,4,5,6,7,8]}
    """
    i = 1
    generic_words = {}
    for w in lengths:
        generic_words[w] = range(i, i + lengths[w])
        i += lengths[w]
    return generic_words

def translate_generic_expression(generic_words, expression):
    """
    Translates an expression involving words into its generic word
    Three special symbols may be used in an expression, -, ~ and ^,
    playing respectively the role of the complement operator, the
    reversal operator and the hat operator. They have to be placed
    in front of the word on which it applies.

    EXAMPLES:
        sage: translate_generic_expression({'u' : [1,2,3,4]}, 'uu')
        [1,2,3,4,1,2,3,4]
        sage: translate_generic_expression({'u' : [1,2,3,4]}, '-u~u^u')
        [-1,-2,-3,-4,4,3,2,1,-4,-3,-2,-1]
    """
    i = 0
    generic_expression = []
    while i < len(expression):
        if expression[i] == '~':
            reversal = copy(generic_words[expression[i+1]])
            reversal.reverse()
            generic_expression += reversal
            i += 2
        elif expression[i] == '-':
            generic_expression += map(lambda n:-n,\
                                      generic_words[expression[i+1]])
            i += 2
        elif expression[i] == '^':
            reversal = copy(generic_words[expression[i+1]])
            reversal.reverse()
            generic_expression += map(lambda n:-n, reversal)
            i += 2            
        else:
            generic_expression += generic_words[expression[i]]
            i += 1
    return generic_expression

def reduce_solution(solution):
    """
    Transforms the solution by reducing the numbers appearing
    For instance, if numbers 1, 3, 4 and 6 appears they are replaced
    respectively by 1, 2, 3 and 4.

    EXAMPLES:
        sage: reduce_solution({'u' : [1,2,2,1], 'v' : [1,5,1]})
        {'u' : [1,2,2,1], 'v' : [1,3,1]}
    """

    numbers = sorted(list(set(map(lambda x:abs(x),\
              reduce(lambda s, t: s + t, solution.values())))))
    reduced = {}
    for w in solution:
        reduced[w] = map(lambda n:(-1) ** int(n < 0) *\
                     (numbers.index(abs(n)) + 1), solution[w])
    return reduced

def overlap(lengths, equations, delays):
    """
    Returns the general solution of the given equation on words

    EXAMPLES::

        sage: overlap({'u' : 4}, ['u', '~u'], [0])
        {'u' : [1,2,2,1]}
        sage: overlap({'u' : 4, 'v' : 4}, ['u', '~u', 'v'], [0, 2])
        {'u' : [1,2,2,1], 'v' : [2,1,3,4]}
    """
    generic_words = get_generic_words(lengths)
    generic_expressions = map(lambda e:translate_generic_expression\
                              (generic_words, e), equations)
    partition = get_generic_partition(generic_expressions, delays)
    if not is_realizable(partition):
        return None
    else:
        representants = get_representants(partition)
        solution = generic_words.copy()
        for gw in solution:
            solution[gw] = map(lambda w:representants[w], generic_words[gw])
        return reduce_solution(solution)

#---------------------------#
# Enumerating all solutions #
#---------------------------#

def overlap_dps_distinct_lengths(length1, length2, length3, delay):
    return overlap({'A' : length1, 'B' : length2, 'X' : length3,\
                    'Y' : length1 + length2 - length3},\
                    ['AB-A-BAB-A-B', '~A~B', 'XY-X-Y', '~X~Y'],\
                    [0, delay, delay])

def get_all_distinct_lengths_dps_particular_solutions(l, m, n):
    solutions = []
    for length1 in range(0, l):
        for length2 in range(0, m):
            for length3 in set(range(0, length1 + length2 + 1)):
                for delay in range(1, length1):
                    solution = overlap_dps_distinct_lengths\
                               (length1, length2, length3, delay)
                    if solution != None:
                        boundary_general = solution['A'] + solution['B']\
                                           + map(lambda n:-n, solution['A'])\
                                           + map(lambda n:-n, solution['B'])
                        if len(set(map(lambda n:abs(n), boundary_general)))\
                           == 2 and seems_prime(boundary_general):
                            solutions.append(Freeman(map(lambda n:PHI(n),\
                                                         boundary_general)))
    return solutions

def get_all_distinct_lengths_dps(l, m, n):
    return filter(lambda p:p[1:].is_simple(),\
           get_all_distinct_lengths_dps_particular_solutions(l, m, n))

def double_square_iterator(max_length=None,\
                           max_iteration=None,\
                           merge_isometric=True):
    i = 0
    p = 2
    double_squares = set([])
    if max_length is None:
        max_length = Infinity
    if max_iteration is None:
        max_iteration = Infinity
    while p <= max_length / 2:
        for a in range(1, p - 1):
            b = p - a
            for d in range(1, min(a - 1, p / 2)):
                for x in range(a + 1 - d, p - 1):
                    y = p - x
                    solution = overlap_dps_distinct_lengths(a, b, x, d)
                    if solution:
                        boundary_general = solution['A'] + solution['B'] + \
                        map(BAR, solution['A']) + map(BAR, solution['B'])
                        if len(set(map(lambda n:abs(n), boundary_general)))\
                           == 2 and seems_prime(boundary_general):
                            word = min(Freeman(map(lambda n:PHI(n),\
                                   boundary_general)).conjugates())
                            if word.is_simple() and word not in double_squares:
                                yield word
                                inv_word = min(F_HAT(word).conjugates())
                                if merge_isometric:
                                    iso_tiles = map(lambda p:\
                                                    min(p.conjugates()),\
                                                    isometric_paths(word))
                                    iso_tiles += map(lambda p:\
                                                     min(p.conjugates()),\
                                                     isometric_paths(inv_word))
                                else:
                                    iso_tiles = [word, inv_word]
                                i += 1
                                if i >= max_iteration:
                                    raise StopIteration
                                double_squares |= set(iso_tiles)
        p += 2

#------------------#
# Useful functions #
#------------------#

def seems_prime(boundary_general):
    i = 0
    prime = True
    while prime and i < len(boundary_general):
        prime = (i % 2 == 0 and boundary_general[i] in [-1, 1])\
                or (i % 2 == 1 and boundary_general[i] in [-2, 2])
        i += 1
    return prime

def plot_tiles(tiles):
    G = [p.plot(pathoptions=dict(rgbcolor='red',thickness=1),\
        endarrow=False,startpoint=False) for p in tiles]
    cols = int(sqrt(len(tiles)))
    rows = len(tiles) / cols + 1
    g = graphics_array(G, rows, cols)
    g.show(figsize=[15, 15 * rows/cols])

def latex_table(tiles, figure_width=2, num_cols=8):
    perimeter_to_tiles = {}
    for tile in tiles:
        if len(tile) not in perimeter_to_tiles:
            perimeter_to_tiles[len(tile)] = [tile]
        else:
            perimeter_to_tiles[len(tile)].append(tile)
    s = '\\begin{supertabular}{|c|l|}\n'
    s += '  \hline\n'
    s += "  P\\'erim\`etre & Tuiles \\\\\n"
    s += '  \hline\n'
    perimeters = sorted(perimeter_to_tiles.keys())
    for p in perimeters:
        s += '  ' + str(p) + ' & \\begin{tikzpicture}\n'
        s += '          \\matrix [ampersand replacement=\&,\n'
        s += '                    every cell/.style={anchor=base},\n'
        s += '                    column sep=5mm,row sep=7mm]\n'
        s += '          {\n'
        i = 0
        for t in perimeter_to_tiles[p]:
            figure_height = figure_width
            h = t.height()
            w = t.width()
            step = 1.0 * figure_width / max(h, w)
            x = min(map(lambda p: p[0], t.points())) * step
            y = min(map(lambda p: p[1], t.points())) * step
            s += '           \\filldraw[draw=black, fill=black!20,\
                                        x=%scm, y=%scm, xshift=%scm,\
                                        yshift=%scm] '\
                 % (step, step, -x + (figure_width - w * step) / 2.0, -y +\
                   (figure_height - h * step) / 2.0)
            s += t.tikz_trajectory() + '; '
            if i % num_cols == num_cols - 1:
                s += '\\\\\n'
            elif i < len(perimeter_to_tiles[p]) - 1:
                s += '\\&\n'
            else:
                s += '\\\\\n'
            i += 1
        s += '          };\n'
        s += '       \\end{tikzpicture}\n'
        s += '       \\\\\n'
        s += '       \hline\n'
    s += '\\end{supertabular}'
    from sage.misc.latex import LatexExpr
    return LatexExpr(s)
        
def isometric_paths(path):
    return [path, F_RHO(path), F_RHO(F_RHO(path)), F_RHO(F_RHO(F_RHO(path))),\
            F_SIGMA(path), F_RHO(F_SIGMA(path)), F_RHO(F_RHO(F_SIGMA(path))),\
            F_RHO(F_RHO(F_RHO(F_SIGMA(path))))]
