#!/usr/bin/env python
# -*- coding: UTF-8 -*-
"""
`extract_doc.py` is a script to extract documentation and
function signatures from Matlab functions and convert them to python-
like input for the Sphinx autodoc system, to generate HTML pages for
API documentation.

Modification for this verion: 

Files are stores separately rather than concatednated to standard out

TODO: there is a bug for functions that have a trailing comma in their
declaration; workaround for now is to manually forbid this as a style
convention.

"""

from __future__ import absolute_import
from __future__ import with_statement
from __future__ import division
from __future__ import nested_scopes
from __future__ import generators
from __future__ import unicode_literals
from __future__ import print_function

import os
import sys
import errno

def ensure_dir(dirname):
    """
    Ensure that a named directory exists; if it does not, attempt to create it.
    http://stackoverflow.com/questions/944536/efficient-way-of-creating-recursive-paths-python
    """
    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            pass

def string(o):
    '''
    try to decode a utf-8 string, might need str or unicode
    depending on python version
    '''
    try:
        return str(o).decode('utf-8')
    except:
        return unicode(o).decode('utf-8')
        
if __name__=='__main__':

    # Optionally specify intput and output directories as command line
    # arguments. For the moment, used to manually handle "submodules"
    indir = '../NFCP/'
    if len(sys.argv)>1:
        indir = sys.argv[1]

    outdir = './_auto/'
    if len(sys.argv)>2:
        outdir = sys.argv[2]
        ensure_dir(outdir)
        # Make sure the target dir is viewed like a module to Python
        outname = '__init__.py'
        with open(outdir+outname,'w') as outfile:
            print('''#!/usr/bin/env python
            # -*- coding: UTF-8 -*-
            ''')
    else:
        ensure_dir(outdir)

    files = [f for f in os.listdir(indir) if f[-2:].lower()=='.m']

    # shadow the print function
    oldprint = print
    
    for f in files:
        # remove .m and add .py extension for putput name
        outname = '.'.join(f.split('.')[:-1])+'.py'
        oldprint('writing',outname)
        with open(outdir+outname,'w') as outfile:
            # shadow print: writes lines to file as well
            # as showing them to standard out. 
            def print(*args):
                line = u' '.join([string(a) for a in args])
                oldprint(line)
                outfile.write(line.encode('utf-8'))
                outfile.write('\n')
            # Declare python source type
            print('#!/usr/bin/env python')
            print('# -*- coding: UTF-8 -*-')
            with open(indir+f) as matlabfile:
                # get ALL lines of tile
                # only works w ASCII
                # lines = matlabfile.read().split(u',')
                lines = [line.decode('utf-8')[:-1] for line in matlabfile.readlines()]
                '''
                Goal:
                Detect all functions
                extract signature
                detect "docstrings", extract.
                '''
                merged_lines       = []
                this_line          = u''
                in_function_header = 0
                # pre-processing pass: remove line continuations that appear
                # in the middle of function declarations (they cause trouble)
                for line in lines:
                    if line.strip()[-3:]!='...':
                        # if no line continuation, emit
                        merged_lines.append(this_line+line.rstrip())
                        this_line = u''
                        in_function_header = 0
                    else:
                        # special edge case: handle line-breaks in the
                        # function declaration itself
                        # if line continuation and in function header, join
                        tokens = line.rstrip().split()                    
                        first = tokens[0]
                        if first==u'function' or in_function_header:
                            in_function_header = 1
                            this_line+=line.rstrip()[:-3]
                # Begin to output source code
                #print(u'# extracted from',f)
                in_comment = 0
                comment    = ''
                in_source  = 0
                source     = ''
                skip_this_source_line = 0
                INDENT     = u''
                for line in merged_lines:
                    # break apart line on whitespace
                    tokens = line.rstrip().split()
                    if len(tokens)>0:
                        first = tokens[0]
                        if not in_comment and first==u'function':
                            INDENT = u'    '
                            tokens = tokens[1:]
                            if any([u'=' in t for t in tokens]):
                                # Need to remove matlab return signature
                                while len(tokens)>0 and not u'=' in tokens[0]:
                                    tokens=tokens[1:]
                                s = tokens[0].strip()
                                if s==u'=':
                                    tokens = tokens[1:]
                                else:
                                    tokens[0] = s.split(u'=')[-1]
                            if in_source:
                                # we've hit a sub-function declaration
                                # emit previous source as quoted comment
                                # (python doc generator will choke on matlab syntax)
                                # print(u"    # begin function source")
                                print(INDENT+u"r'''#STARTCODE")
                                print(source.encode('utf-8'))
                                print(INDENT+u"'''#STOPCODE")
                                # print(u"    # end function source")
                            print(u'def',u' '.join(tokens)+':')
                            in_source = 1
                            source = ''
                            skip_this_source_line = 1
                        if first==u'%{':
                            # Begin block comment (header?)
                            in_comment = 1
                            line = u' '.join(tokens[1:])
                            # print(u"    # begin block comment")
                            print(INDENT+u"r'''")
                        if first==u'%}':
                            # End block comment (header?)
                            in_comment=0
                            print(comment.encode('utf-8'))
                            print(INDENT+u"'''")
                            # print(u"    # end block comment")
                            print(INDENT+u"pass#SKIPME")
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
                print(INDENT+u"'''#STARTCODE")
                print(source.encode('utf-8'))
                print(INDENT+u"'''#STOPCODE")
