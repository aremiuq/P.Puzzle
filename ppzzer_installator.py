
from distutils.sysconfig import get_python_lib
from setuptools import find_packages, setup
from distutils.core import setup
import sys
import os

overlay_warning = False
if "install" in sys.argv:
    lib_paths = [get_python_lib()]
    if lib_paths[0].startswith("/usr/lib/"):
        # We have to try also with an explicit prefix of /usr/local in order to
        # catch Debian's custom user site-packages directory.
        lib_paths.append(get_python_lib(prefix="/usr/local"))
    for lib_path in lib_paths:
        existing_path = os.path.abspath(os.path.join(lib_path, "ppzzer"))
        if os.path.exists(existing_path):
            # We note the need for the warning here, but present it after the
            # command is run, so it's more likely to be seen.
            overlay_warning = True
            break

setup(name="Protein Puzzler",
      version="1.0",
      author="Aleix Matabacas & Austin McKitrick",
      author_email="aleix.matabacas01@estudiant.upf.edu",
      description="Develops the quaternary structure of a protein from their interaction pairs",
      license="BSD",
      packages=find_packages(),
      include_package_data=True,
      scripts=["ppzzer.py"],
      install_requires=["biopython"],
      zip_safe=False,
     )

if overlay_warning:
    sys.stderr.write("""
========
WARNING!
========
You have just installed ppzzr over top of an existing
installation, without removing it first. Because of this,
your install may now include extraneous files from a
previous version that have since been removed from
ppzzr. This is known to cause a variety of problems. You
should manually remove the
%(existing_path)s
directory and re-install ppzzr.
""" % {"existing_path": existing_path})
