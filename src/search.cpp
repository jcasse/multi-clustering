#include <ctime>                // provides: time, clock, CLOCKS_PER_SEC
#include "Multiclustering.h"

using namespace std;

namespace rlair_multi_clustering
{
  bool Multiclustering::add_cluster(int way)
  // Library facilities used: none
  // currently, this function tries each unit in order, but this is not
  // deterministic because it depends on the order in which the data is.
  // to make it deterministic, the function could be changed to compute
  // the resulting cost of removing each unit and selecting the unit
  // decreses the cost the most and repeat. this is very expensive, but maybe
  // the relative order of cost remains after removing units and the cost does
  // not need to be recalculated each time before removing the next unit.
  {
    // select cluster to split
    int cluster = split_cluster(way);
    if(cluster == -1) return false; // clusters are perfect

    // initialize
    cluster_t new_cluster;
    int & values = data->values;
    clustering_t & clustering = clusterings[way];
    cluster_t & cluster_struct = clustering[cluster];
    int units = int(cluster_struct.size());

    // get units signatures and block signature
    int blocks = blocking_size(way);
    vector<vector<counts_t> > units_signatures
      (units, vector<counts_t>(blocks, counts_t(values)));
    for(int i = 0; i != units; ++i)
      get_unit_signature(units_signatures[i], way, cluster, i);
    vector<counts_t> cluster_signature(blocks, counts_t(values));
    for(int b = 0; b != blocks; ++b)
      for(int v = 0; v != values; ++v)
        for(int u = 0; u != units; ++u)
          cluster_signature[b][v] += units_signatures[u][b][v];

    // initial average cluster cost
    double average_cluster_cost = cluster_cost(cluster_signature) / units;

    // move units (from the crossassociation paper)
    int current_units = units;
    for(int i = units - 1; i >= 0; --i)
    {
      // update block counts minus unit counts
      vector<counts_t> block_counts = cluster_signature;
      for(int b = 0; b != blocks; ++b)
        for(int v = 0; v != values; ++v)
          block_counts[b][v] -= units_signatures[i][b][v];
      double cost = cluster_cost(block_counts) / (current_units - 1);
      
      // move unit into new cluster and remove from current cluster
      if(cost < average_cluster_cost)
      {
        new_cluster.push_back(cluster_struct[i]);
        cluster_struct.erase(cluster_struct.begin() + i);
        average_cluster_cost = cost;
        cluster_signature = block_counts;
        --current_units;
      }
    }

    // check status
    if(new_cluster.empty()) return false;

    // add new cluster
    clustering.push_back(new_cluster);

    // erase old cluster if empty
    if(clustering[cluster].empty())
      clustering.erase(clustering.begin() + cluster);

    return true;
  }

  int Multiclustering::split_cluster(int way)
  // Library facilities used: none
  {
    int index = -1;
    double highest_cost = 0;
    int clusters = int(clusterings[way].size());
    for(int cluster = 0; cluster < clusters; cluster++)
    {
      int units = int(clusterings[way][cluster].size());
      double average_cluster_cost = cluster_cost(way, cluster) / units;
      if(average_cluster_cost > highest_cost)
      {
        highest_cost = average_cluster_cost;
        index = cluster;
      }
    }
    return index;
  }

  void Multiclustering::get_unit_signature
  (vector<counts_t> & signature, int way, int cluster, int unit_index)
  // Library facilities used: none
  {
    int & values = data->values;
    vector<int> dimensions = data->dimensions();
    Indexer indexer = blocking_indexer(way, cluster);
    while(!indexer.end())
    {
      // compute value counts for way unit in block
      vector<int> block_counts(values);
      multicluster_t block = get_block(indexer.get_tuple());
      Indexer unit_indexer = block_indexer(block, way, unit_index);
      while(!unit_indexer.end())
      {
        vector<int> tuple = unit_indexer.get_tuple();
        ++block_counts[data->matrix[hyper_index(tuple, dimensions)]];
        unit_indexer.forward();
      }

      // add to unit counts
      int index = indexer.get_sub_index();
      signature[index] = block_counts;

      // index
      indexer.forward();
    }
  }

  double Multiclustering::cluster_cost(const vector<counts_t> & block_counts)
  // Library facilities used: none
  {
    double cost = 0;
    for(size_t b = 0; b != block_counts.size(); ++b)
      cost += hoffman_coding(block_counts[b]);
    return cost;
  }

  double Multiclustering::cluster_cost(int way, int cluster)
  // Library facilities used: none
  {
    double cost = 0;
    Indexer indexer = blocking_indexer(way, cluster);
    while(!indexer.end())
    {
      cost += block_cost(indexer.get_tuple());
      indexer.forward();
    }
    return cost;
  }
}
