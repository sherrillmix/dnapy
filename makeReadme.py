import pkg_resources
import subprocess
import dnapy
import re

entries=pkg_resources.get_entry_map('dnapy',None)['console_scripts'].keys()

helps=[subprocess.check_output(x+' --help;exit 0',stderr=subprocess.STDOUT,shell=True) for x in entries]
print(entries)
#print(helps)
inserts=['\n'.join([head,"~~~~","","::",'  ']+['  '+x for x in help.decode('utf8').split('\n')]) for head,help in zip(entries,helps)]

insert='\n'.join(inserts)

with open("preREADME.rst","r") as pre:
    with open("README.rst","w") as read:
        read.write(pre.read().replace("INSERT_USAGE_HERE",insert))


