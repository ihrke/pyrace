import os, glob
fnames=glob.glob('notebooks/*')

template="http://nbviewer.ipython.org/urls/raw.github.com/ihrke/pyrace/notebook_tests/notebooks/{notebook}"

for fname in fnames:
	fname=os.path.basename(fname)
	print '* [%s](%s)'%(os.path.splitext(fname)[0],template.format(notebook=fname))
