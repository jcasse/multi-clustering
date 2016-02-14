// FILE: HyperMatrix.h
// TEMPLATE CLASS PROVIDED: HyperMatrix (part of the namespace rlair_multi_clustering)

#ifndef RLAIR_MULTI_CLUSTERING_HYPERMATRIX
#define RLAIR_MULTI_CLUSTERING_HYPERMATRIX

#include <fstream>                  // provides: ostream
#include <string>                   // provides: string
#include <vector>                   // provides: vector
#include "Indexer.h"

namespace rlair_multi_clustering
{
  template<class T>
  class HyperMatrix
  {
  public:
    friend class Data;

    // CONSTRUCTORS and DESTRUCTORS
    HyperMatrix() {}

    HyperMatrix(const HyperMatrix & source)
    : data(source.data), dimensions(source.dimensions) {}

    HyperMatrix(const std::vector<int> & dimensions)
    : dimensions(dimensions)
    // Library facilities used: none
    {
      size_t size = 1;
      for(size_t i = 0; i != dimensions.size(); ++i) size *= dimensions[i];
      data = std::vector<T>(size);
    }

    // MODIFICATION MEMBER FUNCTIONS
    HyperMatrix & operator=(const HyperMatrix & source)
    // Precondition: none
    // Postcondition: *this == source
    // Library facilities used: none
    {
      data = source.data;
      dimensions = source.dimensions;
      return *this;
    }

    T & operator[](int index)
    // Precondition: index is within range
    // Postcondition: Return value is the element located by flat index
    // Library facilities used: none
    {
      assert(index < static_cast<int>(data.size()));
      return data[index];
    }

    T operator[](int index) const
    // Precondition: index is within range
    // Postcondition: Return value is the element located by flat index
    // Library facilities used: none
    {
      assert(index < static_cast<int>(data.size()));
      return data[index];
    }

    // CONSTANT MEMBER FUNCTIONS
    size_t size() const {return data.size();}
    // Precondition: none
    // Postcondition: The return value is the number of elements in matrix

    bool operator==(const HyperMatrix & rhs) const
    // Precondition: none
    // Postcondition: Return value is *this == rhs
    {return data == rhs.data && dimensions == rhs.dimensions;}

    bool operator!=(const HyperMatrix & rhs) const
    // Precondition: none
    // Postcondition: Return value is *this != rhs
    {return !operator==(rhs);}

    void print_2D_slice
    (const std::vector<int> & dimension, const std::string & file) const
    // Precondition: file is the name of output file
    // Postcondition: matrix.data is printed to file
    // Library facilities used: assert
    {
      std::ofstream out(file.c_str());
      print_2D_slice(dimension, out);
      out.close();
    }

    void print_2D_slice
    (const std::vector<int> & dimension, std::ostream & out) const
    // Precondition: out is an open output stream
    // Postcondition: matrix.data is printed to out
    // Library facilities used: assert
    {
      int ways = static_cast<int>(dimensions.size());

      assert(int(dimension.size()) == ways);

      // get row and column ways
      std::vector<int> plane;
      for(int way = 0; way != ways; ++way)
        if(dimension[way] == -1) plane.push_back(way);
      assert(int(plane.size()) == 2);

      const int ROW = plane[0];
      const int COL = plane[1];

      assert(ROW < COL);

      // set indexer tuple
      std::vector<int> tuple(ways);
      for(int way = 0; way != ways; ++way)
        if(dimension[way] == -1) tuple[way] = 0;

      // build hyper-cube indexer
      std::vector<bool> mask(ways, true); mask[ROW] = mask[COL] = false;
      Indexer * indexer = new Indexer(*this, tuple, mask);

      // print
      int ccount = 0;
      while(!indexer->end())
      {
        out << operator[](indexer->index());
        if(++ccount % dimensions[COL] == 0) out << " " << std::endl;
        indexer->forward();
      }
    }

    std::vector<std::vector<int> > indexes() const
    // Library facilities used: none
    {
      int ways = int(dimensions.size());
      std::vector<std::vector<int> > temp(ways);
      for(int way = 0; way != ways; ++way)
        for(int index = 0; index != dimensions[way]; ++index)
          temp[way].push_back(index);
      return temp;
    }

    // MEMBER VARIABLES
    std::vector<T> data;            // complete data
    std::vector<int> dimensions;    // size of each way
  };
}

//#include "HyperMatrix.cpp"

#endif
