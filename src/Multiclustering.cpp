#include <ctime>                // provides: time, clock, CLOCKS_PER_SEC
#include <iomanip>              // provides: setw
#include <cmath>                // provides: log
#include <algorithm>            // provides: sort
#include "Indexer.h"
#include "Multiclustering.h"

using namespace std;

namespace rlair_multi_clustering
{
  Multiclustering::Multiclustering() : data(NULL), lout(NULL) {}

  Multiclustering::Multiclustering(Data * data, ostream * lout)
  : data(data), lout(lout)
  // Library facilities used: none
  {initialize();}
  
  Multiclustering::Multiclustering
  (Data * data, std::ostream * lout, std::vector<int> & clusters)
  : data(data), lout(lout)
  // Library facilities used: none
  {initialize(clusters);}

  Multiclustering::Multiclustering(const Multiclustering & source)
  : data(source.data), lout(source.lout)
  // Library facilities used: none
  {copy(source);}

  void Multiclustering::copy(const Multiclustering & source)
  // Library facilities used: none
  {
    data = source.data;
    clusterings = source.clusterings;
    lout = source.lout;
  }

  void Multiclustering::initialize()
  // Library facilities used: none
  {
    // make multi-clustering of single cluster per mode
    vector<int> clusters(data->ways(), 1);
    initialize(clusters);
  }
  
  void Multiclustering::initialize(std::vector<int> & clusters)
  // Library facilities used: none
  {make_multiclustering(clusters);}

  Multiclustering & Multiclustering::operator=(const Multiclustering & source)
  // Library facilities used: none
  {copy(source); return *this;}

  void Multiclustering::make_multiclustering(vector<int> & clusters)
  // Library facilities used: vector
  {
    int ways = data->ways();
    vector<int> cluster(ways);
    clusterings = vector<clustering_t>(ways);
    for(int way = 0; way != ways; ++way)
    {
      clusterings[way] = clustering_t(clusters[way]);
      for(int i = 0; i != data->matrix.dimensions[way]; ++i)
      {
        clusterings[way][cluster[way]].push_back(i);
        cluster[way] = (cluster[way] + 1) % clusters[way];
      }
    }
  }

  int Multiclustering::blocking_size() const
  // Library facilities used: none
  {
    int size = 1;
    for(int i = 0; i != data->ways(); ++i)
      size *= int(clusterings[i].size());
    return size;
  }

  int Multiclustering::blocking_size(int way) const
  // Library facilities used: none
  {
    int size = 1;
    for(int i = 0; i != data->ways(); ++i)
      if(i != way) size *= int(clusterings[i].size());
    return size;
  }


  Indexer Multiclustering::blocking_indexer(int way, int cluster)
  // Library facilities used: assert
  {
    int ways = data->ways();
    assert(way < ways);
    assert(cluster < int(clusterings[way].size()));

    // create indexes
    Indexer::indexes_t indexes(ways);
    for(int i = 0; i != ways; ++i)
      for(int j = 0; j != int(clusterings[i].size()); ++j)
        indexes[i].push_back(j);

    // set tuple
    vector<int> tuple(ways);
    tuple[way] = cluster;

    // set mask for way
    Indexer::mask_t mask(ways);
    mask[way] = true;

    return Indexer(indexes, tuple, mask);
  }

  Indexer Multiclustering::block_indexer
  (multicluster_t & multicluster, int way, int unit_index)
  // Library facilities used: assert
  {
    int ways = data->ways();
    assert(int(multicluster.size()) == ways);
    assert(way < ways);
    assert(unit_index < int(multicluster[way].size()));

    // create indexes
    Indexer::indexes_t indexes(ways);
    for(int i = 0; i != ways; ++i) indexes[i] = multicluster[i];

    // set tuple
    Indexer::tuple_t tuple(ways);
    tuple[way] = unit_index;

    // set mask for way
    Indexer::mask_t mask(ways);
    mask[way] = true;

    return Indexer(indexes, tuple, mask);
  }

  multicluster_t Multiclustering::get_block(const Indexer::tuple_t & tuple)
  // Library facilities used: assert
  {
    int ways = data->ways();
    multicluster_t multicluster(ways);
    for(int way = 0; way != ways; ++way)
      multicluster[way] = clusterings[way][tuple[way]];
    return multicluster;
  }

  Indexer::dimensions_t Multiclustering::blocking_dimensions() const
  // Library facilities used: none
  {
    Indexer::dimensions_t dimensions(data->ways());
    for(int way = 0; way != data->ways(); ++way)
      dimensions[way] = int(clusterings[way].size());
    return dimensions;
  }

  int Multiclustering::block_size(const Indexer::tuple_t & tuple) const
  // Library facilities used: none
  {
    int ways = data->ways();
    assert(int(tuple.size()) == ways);
    int size = 1;
    for(int way = 0; way != ways; ++way)
      size *= int(clusterings[way][tuple[way]].size());
    return size;
  }
}
