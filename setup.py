from setuptools import setup, find_packages

setup(name="mifeos",
      version="1.0.0",
      description="minkowski functional estimator",
      license="GPLv3",
      url="https://github.com/gioacchinowang/mifeos",
      maintainer="Jiaxin Wang",
      maintainer_email="jiaxin.wang@sjtu.edu.cn",
      packages=find_packages(),
      include_package_data=True,
      platforms="any",
      python_requires='>=3.5',
      install_requires=['numpy'],
      zip_safe=False,
      classifiers=["Development Status :: 5 - Production/Stable",
                   "Topic :: Utilities",
                   "License :: OSI Approved :: GNU General Public License v3 "
                   "or later (GPLv3+)"],)
