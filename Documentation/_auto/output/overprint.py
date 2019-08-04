#!/usr/bin/env python
# -*- coding: UTF-8 -*-
def overprint(msg,ntoclear):
    r'''

    Erase current line then print message

    Parameters
    ----------
    msg : `string`
        Message to print
    ntoclear : `int`
        Number of characters to clear

    Returns
    -------
    ntoclear : `int`
        Length of `msg`; number of characters to clear on next iteration

    '''
    pass#SKIPME
    '''#STARTCODE
	fprintf(repmat('\b', 1, ntoclear));
	fprintf(msg);
	ntoclear = numel(msg);

    '''#STOPCODE
