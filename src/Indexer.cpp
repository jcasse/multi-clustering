// FILE: Indexer.cpp (part of namespace rlair_multi_clustering)
// TEMPLATE CLASS implemented: Hyperator (see Indexer.h for documentation)

#include <cassert>                  // provides: assert
#include "Indexer.h"
#include "HyperMatrix.h"

using namespace std;

namespace rlair_multi_clustering
{
  Indexer::Indexer() {}

  Indexer::Indexer(indexes_t & indexes, mask_t & mask)
    : indexes(indexes), mask(mask)
  // Library facilities used: assert
  {
    assert(indexes.size() == mask.size());

    // set dimensions
    size_t ways = indexes.size();
    dimensions = dimensions_t(ways);
    for(size_t way = 0; way != ways; ++way)
      dimensions[way] = static_cast<int>(indexes[way].size());

    // initialize tuple
    tuple = tuple_t(ways);
  }

  Indexer::Indexer(indexes_t & indexes, tuple_t & tuple, mask_t & mask)
    : indexes(indexes), tuple(tuple), mask(mask)
  // Library facilities used: assert
  {
    assert(indexes.size() == tuple.size());
    assert(tuple.size() == mask.size());

    // set dimensions
    size_t ways = indexes.size();
    dimensions = dimensions_t(ways);
    for(size_t way = 0; way != ways; ++way)
      dimensions[way] = static_cast<int>(indexes[way].size());
  }

  Indexer::Indexer(dimensions_t & dimensions, mask_t & mask)
    : dimensions(dimensions), mask(mask)
  // Library facilities used: assert
  {
    assert(dimensions.size() == mask.size());

    // build default indexes
    indexes = default_indexes();

    // initialize tuple
    tuple = tuple_t(dimensions.size());
  }

  Indexer::Indexer
  (const dimensions_t & dimensions, indexes_t & indexes, mask_t & mask)
    : dimensions(dimensions), indexes(indexes), mask(mask)
  // Library facilities used: none
  {
    assert(dimensions.size() == indexes.size());
    assert(indexes.size() == mask.size());

    // initialize tuple
    tuple = tuple_t(dimensions.size());
  }
  
  Indexer::Indexer
  (const HyperMatrix<T> & matrix, const tuple_t & tuple, const mask_t & mask)
    : tuple(tuple), mask(mask)
  {
    dimensions = matrix.dimensions;
    indexes = matrix.indexes();
  }

  Indexer::~Indexer() {}

  void Indexer::reset()
  // Library facilities used: none
  {
    for(size_t way = 0; way != tuple.size(); ++way)
      if(!mask[way]) tuple[way] = 0;
  }

  void Indexer::set(const tuple_t & tuple)
  // Library facilities used: none
  {
    assert(this->tuple.size() == tuple.size());
    this->tuple = tuple;
  }

  void Indexer::forward()
  // Library facilities used: none
  {
    // get dimensionality
    int ways = static_cast<int>(dimensions.size());

    // check if last item
    bool last = true;
    for(int i = 0; i != ways; ++i)
      if(!mask[i]) if(tuple[i] != dimensions[i] - 1) last = false;
    if(last)
    {
      // increment last unmasked index past last item
      int i = ways - 1;
      while(i >=0 && mask[i]) --i;
      assert(i >= 0);
      ++tuple[i];
      return;
    }

    // normal operation
    bool carry = true;
    for(int i = ways - 1; i >= 0; --i)
    {
      if(mask[i]) continue;
      if(!carry) continue;

      // reset index and continue
      if(tuple[i] == dimensions[i] - 1) {tuple[i] = 0;}

      // increment index and stop
      else {++tuple[i]; carry = false;}
    }
  }

  Indexer::indexes_t Indexer::default_indexes()
  // Library facilities used: none
  {
    // get dimensionalty
    size_t ways = dimensions.size();

    // generate full spectrum indexes
    indexes_t indexes(ways);
    for(size_t way = 0; way != ways; ++way)
      for(int i = 0; i != dimensions[way]; ++i)
        indexes[way].push_back(i);

    return indexes;
  }

  int Indexer::index() const
  // Library facilities used: none
  {
    return hyper_index(get_tuple(), dimensions);
  }

  bool Indexer::end() const
  // Library facilities used: assert
  {
    // check last unmasked index
    int i = static_cast<int>(dimensions.size()) - 1;
    while(i >=0 && mask[i]) --i;
    assert(i >= 0);
    return tuple[i] == dimensions[i];
  }

  Indexer::tuple_t Indexer::get_tuple() const
  // Library facilities used: none
  {
    tuple_t indexer_tuple = tuple;
    for(size_t i = 0; i != tuple.size(); ++i)
      indexer_tuple[i] = indexes[i][tuple[i]];
    return indexer_tuple;
  }

  Indexer::tuple_t Indexer::get_sub_tuple() const
  // Library facilities used: none
  {
    tuple_t sub_tuple;
    for(size_t i = 0; i != tuple.size(); ++i)
      if(!mask[i]) sub_tuple.push_back(indexes[i][tuple[i]]);
    return sub_tuple;
  }

  Indexer::dimensions_t Indexer::get_sub_dimensions() const
  // Library facilities used: none
  {
    dimensions_t sub_dimensions;
    for(size_t i = 0; i != dimensions.size(); ++i)
      if(!mask[i]) sub_dimensions.push_back(dimensions[i]);
    return sub_dimensions;
  }

  int Indexer::get_sub_index() const
  // Library facilities used: none
  {
    tuple_t sub_tuple = get_sub_tuple();
    dimensions_t sub_dimensions = get_sub_dimensions();
    return hyper_index(sub_tuple, sub_dimensions);
  }
}
