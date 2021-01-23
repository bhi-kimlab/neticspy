========
Overview
========

.. start-badges

.. list-table::
    :stub-columns: 1

    * - docs
      - |docs|
    * - tests
      - | |travis| |appveyor| |requires|
        | |codecov|
    * - package
      - | |version| |wheel| |supported-versions| |supported-implementations|
        | |commits-since|
.. |docs| image:: https://readthedocs.org/projects/pyNetICS/badge/?style=flat
    :target: https://readthedocs.org/projects/pyNetICS
    :alt: Documentation Status

.. |travis| image:: https://api.travis-ci.org/bhi-kimlab/pyNetICS.svg?branch=master
    :alt: Travis-CI Build Status
    :target: https://travis-ci.org/bhi-kimlab/pyNetICS

.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/github/bhi-kimlab/pyNetICS?branch=master&svg=true
    :alt: AppVeyor Build Status
    :target: https://ci.appveyor.com/project/bhi-kimlab/pyNetICS

.. |requires| image:: https://requires.io/github/bhi-kimlab/pyNetICS/requirements.svg?branch=master
    :alt: Requirements Status
    :target: https://requires.io/github/bhi-kimlab/pyNetICS/requirements/?branch=master

.. |codecov| image:: https://codecov.io/gh/bhi-kimlab/pyNetICS/branch/master/graphs/badge.svg?branch=master
    :alt: Coverage Status
    :target: https://codecov.io/github/bhi-kimlab/pyNetICS

.. |version| image:: https://img.shields.io/pypi/v/neticspy.svg
    :alt: PyPI Package latest release
    :target: https://pypi.org/project/neticspy

.. |wheel| image:: https://img.shields.io/pypi/wheel/neticspy.svg
    :alt: PyPI Wheel
    :target: https://pypi.org/project/neticspy

.. |supported-versions| image:: https://img.shields.io/pypi/pyversions/pynetics.svg
    :alt: Supported versions
    :target: https://pypi.org/project/pynetics

.. |supported-implementations| image:: https://img.shields.io/pypi/implementation/pynetics.svg
    :alt: Supported implementations
    :target: https://pypi.org/project/pynetics

.. |commits-since| image:: https://img.shields.io/github/commits-since/bhi-kimlab/pyNetICS/v0.0.5.svg
    :alt: Commits since latest release
    :target: https://github.com/bhi-kimlab/pyNetICS/compare/v0.0.5...master



.. end-badges

Python implementation of NetICS.

* Free software: MIT license

Installation
============

::

    pip install neticspy

You can also install the in-development version with::

    pip install https://github.com/bhi-kimlab/neticspy/archive/master.zip


Documentation
=============


https://pyNetICS.readthedocs.io/


Development
===========

To run all the tests run::

    tox

Note, to combine the coverage data from all the tox environments run:

.. list-table::
    :widths: 10 90
    :stub-columns: 1

    - - Windows
      - ::

            set PYTEST_ADDOPTS=--cov-append
            tox

    - - Other
      - ::

            PYTEST_ADDOPTS=--cov-append tox
