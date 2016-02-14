// FILE: HyperMatrix.cpp (part of namespace rlair_multi_clustering)
// TEMPLATE CLASS implemented: HyperMatrix (see HyperMatrix.h for documentation)

#include <cassert>                  // provides: assert
#include "Indexer.h"

using namespace std;

namespace rlair_multi_clustering
{
  template <class T>
  HyperMatrix<T>::HyperMatrix() {}

  template <class T>
  HyperMatrix<T>::HyperMatrix(const HyperMatrix & source)
    : data(source.data), dimensions(source.dimensions) {}

  template <class T>
  HyperMatrix<T>::HyperMatrix(const vector<int> & dimensions)
    : dimensions(dimensions)
  // Library facilities used: none
  {
    size_t size = 1;
    for(size_t i = 0; i != dimensions.size(); ++i) size *= dimensions[i];
    data = vector<T>(size);
  }

  template <class T>
  HyperMatrix<T> & HyperMatrix<T>::operator=(const HyperMatrix<T> & source)
  // Library facilities used: none
  {
    data = source.data;
    dimensions = source.dimensions;
    return *this;
  }
  
  template <class T>
  T & HyperMatrix<T>::operator[](int index)
  // Library facilities used: none
  {
    assert(index < static_cast<int>(data.size()));
    return data[index];
  }

  template <class T>
  T HyperMatrix<T>::operator[](int index) const
  // Library facilities used: none
  {
    assert(index < static_cast<int>(data.size()));
    return data[index];
  }

  template <class T>
  size_t HyperMatrix<T>::size() const
  // Library facilities used: none
  {
    return data.size();
  }

  template <class T>
  bool HyperMatrix<T>::operator==(const HyperMatrix<T> & rhs) const
  // Library facilities used: none
  {
    return data == rhs.data && dimensions == rhs.dimensions;
  }

  template <class T>
  bool HyperMatrix<T>::operator!=(const HyperMatrix<T> & rhs) const
  // Library facilities used: none
  {
    return !operator==(rhs);
  }

  template <class T>
  void HyperMatrix<T>::print_2D
  (const std::vector<int> & dimension, const std::string & file) const
  // Library facilities used: assert
  {
    ofstream out(file.c_str());
    print_2D(dimension, out);
    out.close();
  }

  template <class T>
  void HyperMatrix<T>::print_2D
  (const std::vector<int> & dimension, ostream& out) const
  // Library facilities used: assert
  {
    int ways = static_cast<int>(dimensions.size());

    assert(int(dimension.size()) == ways);

    // get row and column ways
    vector<int> plane;
    for(int way = 0; way != ways; ++way)
      if(dimension[way] == -1) plane.push_back(way);
    assert(int(plane.size()) == 2);

    const int ROW = plane[0];
    const int COL = plane[1];

    assert(ROW < COL);

    // set indexer tuple
    vector<int> tuple(ways);
    for(int way = 0; way != ways; ++way)
      if(dimension[way] == -1) tuple[way] = 0;

    // build hyper-cube indexer
    std::vector<bool> mask(ways, true); mask[ROW] = mask[COL] = false;
    Indexer indexer(*this, tuple, mask);

    // print
    int ccount = 0;
    while(!indexer.end())
    {
      out << operator[](indexer.index());
      if(++ccount % dimensions[COL] == 0) out << " " << endl;
      indexer.forward();
    }
  }

  template <class T>
  std::vector<std::vector<int> > HyperMatrix<T>::indexes() const
  // Library facilities used: none
  {
    int ways = int(dimensions.size());
    vector<vector<int> > temp(ways);
    for(int way = 0; way != ways; ++way)
      for(int index = 0; index != dimensions[way]; ++index)
        temp[way].push_back(index);
    return temp;
  }
}
