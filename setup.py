#!/usr/bin/env python
# -*- Mode: python; tab-width: 4; indent-tabs-mode:nil; coding:utf-8 -*-

from setuptools import setup

name = "psntools"
url = "https://www.github.com/ELELAB/psntools"
author = "Valentina Sora, Matteo Tiberti, Elena Papaleo"
author_email = "sora.valentina1@gmail.com"
version = "0.1"
description = "PSN tools"
package_data = {"psntools" : ["config_plot/*"]}
package_dir = {"psntools" : "psntools"}
packages = ["psntools"]
install_requires = \
    ["MDAnalysis", "networkx", "numpy", "pyyaml", "Cython", \
     "matplotlib", "seaborn"]

setup(name = name,
      url = url,
      author = author,
      author_email = author_email,
      version = version,
      description = description,
      package_data = package_data,
      package_dir = package_dir,
      packages = packages,
      install_requires = install_requires)