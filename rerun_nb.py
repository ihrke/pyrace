#!/usr/bin/env python
"""
"""
import io
import os,sys,time
import base64
import re

from collections import defaultdict
from Queue import Empty

try:
    from IPython.kernel import KernelManager
except ImportError:
    from IPython.zmq.blockingkernelmanager import BlockingKernelManager as KernelManager
from IPython.nbformat.current import reads, NotebookNode, writes


def run_cell(shell, iopub, cell):
    # print cell.input
    shell.execute(cell.input)
    # wait for finish, maximum 20s
    shell.get_msg()
    outs = []
    
    while True:
        try:
            msg = iopub.get_msg(timeout=20)
        except Empty:
            break
        msg_type = msg['msg_type']
        if msg_type in ('status', 'pyin'):
            continue
        elif msg_type == 'clear_output':
            outs = []
            continue
        
        content = msg['content']
        # print msg_type, content
        out = NotebookNode(output_type=msg_type)
        
        if msg_type == 'stream':
            out.stream = content['name']
            out.text = content['data']
        elif msg_type in ('display_data', 'pyout'):
            for mime, data in content['data'].iteritems():
                attr = mime.split('/')[-1].lower()
                # this gets most right, but fix svg+html, plain
                attr = attr.replace('+xml', '').replace('plain', 'text')
                setattr(out, attr, data)
            if msg_type == 'pyout':
                #out.prompt_number = content['execution_count']
                #TODO: need to find better workaround
                pass
        elif msg_type == 'pyerr':
            out.ename = content['ename']
            out.evalue = content['evalue']
            out.traceback = content['traceback']
        else:
            print "unhandled iopub msg:", msg_type
        
        outs.append(out)
    return outs
    

def run_notebook(nb):
    km = KernelManager()
    km.start_kernel(extra_arguments=['--pylab=inline'], stderr=open(os.devnull, 'w'))
    try:
        kc = km.client()
        kc.start_channels()
        iopub = kc.iopub_channel
    except AttributeError:
        # IPython 0.13
        kc = km
        kc.start_channels()
        iopub = kc.sub_channel
    shell = kc.shell_channel
    
    # run %pylab inline, because some notebooks assume this
    # even though they shouldn't
    shell.execute("pass")
    shell.get_msg()
    while True:
        try:
            iopub.get_msg(timeout=1)
        except Empty:
            break
    
    successes = 0
    failures = 0
    errors = 0
    prompt_number = 1
    for ws in nb.worksheets:
        for cell in ws.cells:
            if cell.cell_type != 'code':
                continue
            try:
                outs = run_cell(shell, iopub, cell)
            except Exception as e:
                print "failed to run cell:", repr(e)
                print cell.input
                errors += 1
                cell.outputs = [e]
                continue
            
            sys.stdout.write('.')

            cell.outputs = outs
            cell.prompt_number = prompt_number
            if cell.outputs:
                cell.outputs[0]['prompt_number'] = prompt_number
            prompt_number += 1

    print
    print "finished running notebook %s" % nb.metadata.name
    kc.stop_channels()
    km.shutdown_kernel()
    del km


import glob

if __name__ == '__main__':
    fnames=glob.glob('notebooks/*.ipynb')

    for fname in fnames:
        with open(fname) as f:
            print "running %s" % fname
            nb = reads(f.read(), 'json')
        run_notebook(nb)

        strnb=writes(nb, 'json')
        with open(fname, 'w') as f:
            f.write(strnb)
