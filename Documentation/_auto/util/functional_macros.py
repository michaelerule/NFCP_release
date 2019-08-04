#!/usr/bin/env python
# -*- coding: UTF-8 -*-
r'''

Functional programming macros for Matlab. `macros.m` should be called as
a script to load these macros as lambda expressions into the current
workspace.

These macros make it easier to define subroutines from within a script,
without creating a new dedicated file. This allows for more compact and
maintainable scripts.

The main use case is in defining objective functions to be passed to
optimization routines. These functions often need to close over scop
(e.g. large datasets), and are often very short and awkward to maintain
in a separate file.

Matlab's restrictions on syntax lead to some limitaions, and these macros
might be more difficult to read in some use cases.

'''
pass#SKIPME
'''#STARTCODE

'''#STOPCODE
