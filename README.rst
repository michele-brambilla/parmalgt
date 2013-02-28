======================================
parmalgt -- A NSPT simulation program
======================================

Quick start guide.
====================

After you download the source code using the commands
::

  git clone git@bitbucket.org:dirkhesse/parmalgt.git
  git submodule init
  git submodule update

you want to generate and run the configure script,
::

  cd parmalgt
  ./autogen.sh
  ./configure CXXFLAGS="-O3 -fopenmp -DSSE -msse2"

(not we set a bunch of compiler flags in the second command). All this
should succeed if you have

* ``autoconf >= 2.59``
* ``automake >= 1.06``

installed. Next you will want to compile the code. The version without
the gradient flow is compiled using
::

  make LocalQuench

while for the version with Wilson you want to use
::

  make LocalQuenchFlow

Both programs should now be ready. Their behaivours are controlled
with the configuration files ``input`` and ``flow_ipt``
respectively. Both contain comments explaining the options. The
execution of the programs will produce data files containing measured
observables, with the extension ``.bindat``. The data is in binary, if
you call the command
::

  od -t f8 Gp.bindat

you will get an output like this::

                     REAL                     IMAGINARY
  0000000     3.000000000000000e+00    6.476300976980079e-17  <- tree level, Measurement 1
  0000020    -6.025411516049490e-19    3.517633715175191e-18  <- O(g), Meas. 1
  0000040    -1.628038858861708e-01   -8.673617379884035e-19  <- O(g^2), Meas. 1
  0000060     3.166703085779656e-04    4.846514929628224e-05  <- ...
                                                              <- tree level, Meas. 2 ...
