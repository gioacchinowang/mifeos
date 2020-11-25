from setuptools import setup, find_packages

setup(name="mifeos",
      version="0.0.1",
      description="Minkowski Functional Estimator on partial sky",
      license="GPLv3",
      url="https://github.com/gioacchinowang/mifeos",
      packages=find_packages(),
      dependency_links=[],
      python_requires='>=3.5',
      install_requires=['numpy', 'healpy'],
      zip_safe=False,
      classifiers=["Development Status :: 4 - Beta",
                   "Topic :: Utilities",
                   "License :: OSI Approved :: GNU General Public License v3 "
                   "or later (GPLv3+)"],)
