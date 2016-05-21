[//]: # (Markdown: dillinger.io/ shows a nice example of Markdown commands with a viewer.)
[//]: # (Comments in Markdown: http://stackoverflow.com/questions/4823468/comments-in-markdown)
[//]: # (C++ Project Structure: http://hiltmon.com/blog/2013/07/03/a-simple-c-plus-plus-project-structure/)
[//]: # (C++ Library Creation: http://www.adp-gmbh.ch/cpp/gcc/create_lib.html)

# Multi-clustering

This program co-clusters multi-dimensional binary data matrices.
A multi-clustering is a generalization of a co-clustering to any number of
dimensions.
A co-clustering is a clustering of both the rows and columns of a
two-dimensional data matrix.

The algorithm clusters all the dimensions simultaneously.
The main advantage of simultaneous clustering is that information from one
clustering may help in the other clustering and vice versa, creating a sort of
sinergy, resulting in a potentially superior clustering (Inderjit S. Dhillon, Subramanyam Mallela, and Dharmendra S. Modha. 2003. Information-theoretic co-clustering. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining (KDD '03)).

Parts of this algorithm apply the strategy described in Deepayan Chakrabarti, Spiros Papadimitriou, Dharmendra S. Modha, and Christos Faloutsos. 2004. Fully automatic cross-associations. In Proceedings of the tenth ACM SIGKDD international conference on Knowledge discovery and data mining (KDD '04).

### Installation

```sh
$ git clone git@github.com:jcasse/multi-clustering.git
$ cd multi-clustering
$ make
```
### Usage

```sh
$ bin/multi-clustering --help
```

License
----

[//]: # (A short snippet describing the license (MIT, Apache, etc.))

[//]: # (http://choosealicense.com/)

Copyright (C) 2016 Juan Casse

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
