# imat
Interface C++ matrix to access GSL C numerical library

---

A (toy) C++ template matrix library to interact with GSL. Only sort, eigen and
covariance functions were implemented here -- consider them as examples of
calling GSL C numerical library for other linear algebra or matrix related
functions.

Migrated to github in 2020.  The code is fininshed and manged by svn around
2008 (a few adaptations later on) -- old, still functional.

Headers only or template library files:

./imat/

Demo or examples (also as tests in developement):

./demos/

---
Compile and run the test:
$ make
$ ./vector\_demo
$ ./matrix\_demo

