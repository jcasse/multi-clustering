// FILE: Indexer.h
// CLASS PROVIDED: Indexer (part of the namespace rlair_multi_clustering)

#ifndef RLAIR_MULTI_CLUSTERING_INDEXER
#define RLAIR_MULTI_CLUSTERING_INDEXER

//#define NDEBUG

#include <vector>                   // provides: vector
#include <cassert>                  // provides: assert

namespace rlair_multi_clustering
{
  template<class T> class HyperMatrix;

  class Indexer
  {
  public:
    typedef int T;
    typedef std::vector<bool> mask_t;
    typedef std::vector<int> tuple_t;
    typedef std::vector<int> dimensions_t;
    typedef std::vector<std::vector<int> > indexes_t;

    // CONSTRUCTORS and DESTRUCTORS
    Indexer();
    Indexer(indexes_t & indexes, mask_t & mask);
    Indexer(indexes_t & indexes, tuple_t & tuple, mask_t & mask);
    Indexer(dimensions_t & dimensions, mask_t & mask);
    Indexer(const dimensions_t & dimensions, indexes_t & indexes, mask_t & mask);
    Indexer(const HyperMatrix<T> & matrix,
      const tuple_t & tuple, const mask_t & mask);
    ~Indexer();

    // MODIFICATION MEMBER FUNCTIONS
    void reset();
    // Precondition: none
    // Postcondition: tuple is the zero vector
    void set(const tuple_t & tuple);
    // Precondition: tuple is within range
    // Postcondition: this->tuple == tuple
    void forward();
    // Precondition: none
    // Postcondition: tuple is incremented
    // indexes are incremented from last to first (like number system)

    // CONSTANT MEMBER FUNCTIONS
    int index() const;
    // Precondition: none
    // Postcondition: Return value is the flat index
    bool end() const;
    // Precondition: none
    // Postcondition: Return value is true if reached past last item
    tuple_t get_tuple() const;
    // Precondition: none
    // Postcondition: Return value is tuple of indexes
    int get_sub_index() const;
    // Precondition: none
    // Postcondition: Return value is flat index of hyper-plane

  //private:
    dimensions_t dimensions;  // number of items in each way
    indexes_t indexes;        // over which to index
    tuple_t tuple;            // current indexer state
    mask_t mask;              // masked indexes will not be iterated

  private:

    // UTILITY FUNCTIONS
    indexes_t default_indexes();
    // Precondition: none
    // Postcondition: Return value is vector of indexes starting at 0
    tuple_t get_sub_tuple() const;
    // Precondition: none
    // Postcondition: Return value is tuple of indexes excluding masked indexs
    dimensions_t get_sub_dimensions() const;
    // Precondition: none
    // Postcondition: Return value is dimensions excluding masked indexs
  };

  inline int hyper_index
  (const std::vector<int> & tuple, const std::vector<int> & dimensions)
  // Library facilities used: assert
  {
    assert(tuple.size() == dimensions.size());
    int ways = static_cast<int>(dimensions.size());
    int index = 0;
    for(int i = 0; i < ways; ++i)
    {
      int prod = 1;
      for(int j = i + 1; j < ways; ++j) prod *= dimensions[j];
      index += prod * tuple[i];
    }
    return index;
  }
}

#endif
