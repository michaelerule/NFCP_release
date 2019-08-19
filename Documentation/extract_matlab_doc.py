#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
`extract_doc.py` is a script to extract documentation and
function signatures from Matlab functions and convert them to python-
like input for the Sphinx autodoc system, to generate HTML pages for
API documentation.
"""

from __future__ import absolute_import
from __future__ import with_statement
from __future__ import division
from __future__ import nested_scopes
from __future__ import generators
from __future__ import unicode_literals
from __future__ import print_function

if __name__=='__main__':

    import os
    import sys

    indir = '../NFCP/'

    files = [f for f in os.listdir(indir) if f[-2:].lower()=='.m']

    print('#!/usr/bin/env python')
    print('# -*- coding: UTF-8 -*-')

    for f in files:
        with open(indir+f) as text_file:
            # only works w ASCII
            # lines = text_file.read().split(u',')
            lines = [line.decode('utf-8')[:-1] for line in text_file.readlines()]
            '''
            Goal:
            Detect all functions
            extract signature
            detect "docstrings", extract.
            '''
            merged_lines = []
            this_line = u''
            in_function_header=0
            for line in lines:
                if line.strip()[-3:]!='...':
                    # if no line continuation, emit
                    merged_lines.append(this_line+line.rstrip())
                    this_line = u''
                    in_function_header = 0
                else:
                    # if line continuation and in function header, join
                    tokens = line.rstrip().split()                    
                    first = tokens[0]
                    if first==u'function' or in_function_header:
                        in_function_header = 1
                        this_line+=line.rstrip()[:-3]
                    
            print(u'# extracted from',f)
            in_comment = 0
            comment = ''
            in_source = 0
            source = ''
            skip_this_source_line = 0
            for line in merged_lines:
                tokens = line.rstrip().split()
                if len(tokens)>0:
                    first = tokens[0]
                    if first==u'function':
                        tokens = tokens[1:]
                        if any([u'=' in t for t in tokens]):
                            #print('### '+line)
                            # Need to remove return signature
                            while len(tokens)>0 and not u'=' in tokens[0]:
                                tokens=tokens[1:]
                            s = tokens[0].strip()
                            if s==u'=':
                                tokens = tokens[1:]
                            else:
                                tokens[0] = s.split(u'=')[-1]
                        if in_source:
                            # we've probably hit a sub-function declaration
                            # emit previous source as quoted comment
                            # (python doc generator will choke on matlab syntax)
                            print(u"    r'''")
                            print(source.encode('utf-8'))
                            print(u"    '''")
                        print(u'def',u' '.join(tokens)+':')
                        in_source = 1
                        source = ''
                        skip_this_source_line = 1
                    if first==u'%{':
                        in_comment = 1
                        line = u' '.join(tokens[1:])
                        print(u"    r'''")
                    if first==u'%}':
                        in_comment=0
                        print(comment.encode('utf-8'))
                        print(u"    '''")
                        print(u"    pass")
                        comment=''
                        skip_this_source_line=1
                if in_source:
                    if skip_this_source_line:
                        skip_this_source_line = 0
                    # don't include bracketd comments
                    elif not in_comment:
                        source += line+u'\n'
                if in_comment:
                    comment+=line.rstrip()+'\n'
            # reached end of file, should emit any remaining code
            print(u"    '''")
            print(source.encode('utf-8'))
            print(u"    '''")
