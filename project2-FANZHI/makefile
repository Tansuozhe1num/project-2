# Makefile for screensaver
# Copyright (c) 2012 the Massachusetts Institute of Technology
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# FANZHI - Modified makefile for quadtree optimization

# The sources we're building
HEADERS = $(wildcard *.h)
PRODUCT_SOURCES = $(filter-out graphic_stuff.c, $(wildcard *.c))

# What we're building
PRODUCT_OBJECTS = $(PRODUCT_SOURCES:.c=.o)
PRODUCT = screensaver
PROFILE_PRODUCT = $(PRODUCT:%=%.prof) #the product, instrumented for gprof

# What we're building with
# OpenCilk compiler - set OPENCILK_PATH to use custom installation
ifdef OPENCILK_PATH
  CXX = $(OPENCILK_PATH)/bin/clang
else
  CXX = clang
endif
CXXFLAGS = -std=gnu99 -Wall -fopencilk
LDFLAGS = -lrt -lm

include ./cilkutils.mk

# Determine which profile--debug or release--we should build against, and set
# CFLAGS appropriately.

ifeq ($(DEBUG),1)
  # We want debug mode.
  CXXFLAGS += -g -O0 -gdwarf-3
else
  # We want release mode.
  CXXFLAGS += -O3 -DNDEBUG
endif


# By default, make the product.
all:		$(PRODUCT)

# How to build for profiling
prof:		$(PROFILE_PRODUCT)

lint:
	python clint.py *.h *.c


# How to clean up
clean:
	$(RM) $(PRODUCT) $(PROFILE_PRODUCT) *.o *.out


# How to compile a C file
%.o:		%.c $(HEADERS)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) -o $@ -c $<

# How to link the product
$(PRODUCT): LDFLAGS += -lXext -lX11
$(PRODUCT):	$(PRODUCT_OBJECTS) graphic_stuff.o
	$(CXX) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $@ $(PRODUCT_OBJECTS) graphic_stuff.o

# How to build the product, instrumented for profiling
$(PROFILE_PRODUCT): CXXFLAGS += -DPROFILE_BUILD -pg
$(PROFILE_PRODUCT): LDFLAGS += -pg
$(PROFILE_PRODUCT): $(PRODUCT_OBJECTS)
	$(CXX)  $(PRODUCT_OBJECTS) $(LDFLAGS) $(EXTRA_LDFLAGS) -o $(PROFILE_PRODUCT)
