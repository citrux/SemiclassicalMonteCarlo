#!/bin/sh

srcdir="src"
files=$(ls $srcdir)
for f in $files; do
    clang-format-3.6 -i -style=file $srcdir/$f
done
