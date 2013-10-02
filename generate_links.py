import os, glob
fnames=glob.glob('notebooks/*')

template="http://nbviewer.ipython.org/urls/raw.github.com/ihrke/pyrace/notebook_tests/notebooks/{notebook}"

linkstr=""
for fname in fnames:
	fname=os.path.basename(fname)
	linkstr+='* [%s](%s)\n'%(os.path.splitext(fname)[0],template.format(notebook=fname))

with open("README.md.tpl", 'r') as f:
    readme=f.read()

with open("README.md", 'w') as f:
    f.write(readme.format(notebooklinks=linkstr))

print ">> generated new README.md"
print ">> content:"
print "----"
print linkstr
print "----"