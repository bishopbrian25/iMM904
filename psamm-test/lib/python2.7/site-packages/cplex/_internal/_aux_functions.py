# --------------------------------------------------------------------------
# File: _aux_functions.py 
# ---------------------------------------------------------------------------
# Licensed Materials - Property of IBM
# 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
# Copyright IBM Corporation 2008, 2015. All Rights Reserved.
#
# US Government Users Restricted Rights - Use, duplication or
# disclosure restricted by GSA ADP Schedule Contract with
# IBM Corp.
# ------------------------------------------------------------------------
"""


"""


from ..exceptions import CplexError
from .. import six
from ..six.moves import map
from ..six.moves import zip


def validate_arg_lengths(env, caller, arg_list):
    """non-public"""
    arg_lengths = [len(x) for x in arg_list]
    max_length = max(arg_lengths)
    for arg_length in arg_lengths:
        if arg_length != 0 and arg_length != max_length:
            raise CplexError(
                "validate_arg_lengths: Inconsistent arguments to " + caller)
    return max_length


def make_ranges(indices):
    """non-public"""
    ranges = []
    i = 0
    j = 0
    while i < len(indices):
        while j < len(indices) - 1 and indices[j + 1] == indices[j] + 1:
            j += 1
        ranges.append((indices[i], indices[j]))
        i = j + 1
        j = i
    return ranges
    

def apply_freeform_two_args(caller, fn, convert, args):
    """non-public"""
    def con(a):
        if isinstance(a, six.string_types):
            return convert(a)
        else:
            return a
    if len(args) == 2:
        if (isinstance(con(args[0]), six.integer_types) and
            isinstance(con(args[1]), six.integer_types)):
            return fn(con(args[0]), con(args[1]))
        else:
            raise CplexError(
                "apply_freeform_two_args: Wrong argument type to " + caller)
    elif len(args) == 1:
        if isinstance(args[0], (list, tuple)):
            retval = []
            for member in map(fn, *zip(*make_ranges(list(map(con, args[0]))))):
                retval.extend(member)
            return retval
        if isinstance(con(args[0]), six.integer_types):
            return fn(con(args[0]), con(args[0]))[0]
        else:
            raise CplexError(
                "apply_freeform_two_args: Wrong argument type to " + caller)
    elif len(args) == 0:
        return fn(0)
    else:
        raise CplexError(
            "apply_freeform_two_args: Wrong number of arguments to " + caller)


def apply_freeform_one_arg(caller, fn, convert, maxval, args):
    """non-public"""
    def con(a):
        if isinstance(a, six.string_types):
            return convert(a)
        else:
            return a
    if len(args) == 2:
        if (isinstance(con(args[0]), six.integer_types) and
            isinstance(con(args[1]), six.integer_types)):
            return [fn(x) for x in range(con(args[0]), con(args[1]) + 1)]
        else:
            raise CplexError(
                "apply_freeform_one_arg: Wrong argument type to " + caller)
    elif len(args) == 1:
        if isinstance(args[0], (list, tuple)):
            return [fn(x) for x in map(con, args[0])]
        elif isinstance(con(args[0]), six.integer_types):
            return fn(con(args[0]))
        else:
            raise CplexError(
                "apply_freeform_one_arg: Wrong argument type to " + caller)
    elif len(args) == 0:
        return apply_freeform_one_arg(caller, fn, convert, 0,
                                      (list(range(maxval)),))
    else:
        raise CplexError(
            "apply_freeform_one_arg: Wrong number of arguments to " + caller)

def apply_pairs(caller, fn, convert, *args):
    """non-public"""
    def con(a):
        if isinstance(a, six.string_types):
            return convert(a)
        else:
            return a
    if len(args) == 2:
        fn([con(args[0])], [args[1]])
    else:
        a1, a2 = zip(*args[0])
        fn(list(map(con, a1)), list(a2))


def delete_set(caller, fn, convert, max_num, *args):
    """non-public"""
    if len(args) == 0:
        for i in range(max_num):
            fn(0)
    elif len(args) == 1:
        if isinstance(convert(args[0]), six.integer_types):
            fn(convert(args[0]))
        else:
            args = list(map(convert, args[0]))
            args.sort()
            for i, a in enumerate(args):
                fn(convert(a) - i)
    elif len(args) == 2:
        delete_set(caller, fn, convert, max_num,
                   list(range(convert(args[0]), convert(args[1]) + 1)))


class _group:
    """Object to contain constraint groups"""

    def __init__(self, gp):
        """Constructor for the _group object

        gp is a list of tuples of length two (the first entry of which
        is the preference for the group (a float), the second of which
        is a tuple of pairs (type, id), where type is an attribute of
        conflict.constraint_type and id is either an index or a valid
        name for the type).

        Example input: [(1.0, ((2, 0),)), (1.0, ((3, 0), (3, 1)))]
        """
        self._gp = gp

        
def make_group(caller, conv, max_num, c_type, *args):
    """Returns a _group object

    input:
    conv    - a function that will convert names to indices
    max_num - number of existing constraints of a given type
    c_type  - constraint type
    args    - arbitrarily many arguments (see description below)

    If args is empty, every constraint/bound is assigned weight 1.0.

    If args is of length one or more, every constraint/bound is assigned
    a weight equal to the float passed in as the first item.

    If args contains additional items, they determine a subset of
    constraints/bounds to be included.  If one index or name is
    specified, it is the only one that will be included.  If two indices
    or names are specified, all constraints between the first and the
    second, inclusive, will be included.  If a sequence of names or
    indices is passed in, all of their constraints/bounds will be
    included.

    See example usage in _subinterfaces.ConflictInterface.
    """
    if len(args) <= 1:
        cons = list(range(max_num))
    if len(args) == 0:
        weight = 1.0
    else:
        weight = args[0]
    if len(args) == 2:
        weight = args[0]
        if isinstance(conv(args[1]), six.integer_types):
            cons = [conv(args[1])]
        else:
            cons = map(conv, args[1])
    elif len(args) == 3:
        cons = list(range(conv(args[1]), conv(args[2]) + 1))
    return _group([(weight, ((c_type, i),)) for i in cons])
