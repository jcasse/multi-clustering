#include <ctime>                // provides: time, clock, CLOCKS_PER_SEC
#include "Multiclustering.h"

using namespace std;

namespace rlair_multi_clustering
{
  bool Multiclustering::optimize(int way)
  // Library facilities used: assert
  // this function assigns each unit to the cluster in which its encoding cost
  // is the lowest. This function can be parallelized.
  {
    assert(way < data->ways());

    // store units signatures and block frequencies
    int & values = data->values;
    int blocks = blocking_size(way);
    int & units = data->matrix.dimensions[way];
    int clusters = int(clusterings[way].size());
    vector<vector<counts_t> > units_signatures
      (units, vector<counts_t>(blocks, counts_t(values)));
    for(int cluster = 0; cluster != clusters; ++cluster)
      for(int i = 0; i != int(clusterings[way][cluster].size()); ++i)
        get_unit_signature
        (units_signatures[clusterings[way][cluster][i]], way, cluster, i);
    vector<vector<counts_t> > clusters_signatures
      (clusters, vector<counts_t>(blocks, counts_t(values)));
    for(int c = 0; c != clusters; ++c)
      for(int b = 0; b != blocks; ++b)
        for(int i = 0; i != int(clusterings[way][c].size()); ++i)
          for(int v = 0; v != values; ++v)
            clusters_signatures[c][b][v] +=
            units_signatures[clusterings[way][c][i]][b][v];

    bool optimized = false;

    // new clustering (assignments)
    vector<int> new_assignments(data->matrix.dimensions[way], -1);

    // optimize way
    for(int cluster = 0; cluster != clusters; ++cluster)
      for(int i = int(clusterings[way][cluster].size()) - 1; i >= 0; --i)
        if(
          optimize(way, cluster, i,
                   units_signatures, clusters_signatures, new_assignments)
        )
          optimized = true;

    if(!optimized) return false;

    // determine number of remaining clusters
    clusters = 0;
    for(int i = 0; i < (int)new_assignments.size(); i++)
      if(new_assignments[i] >= clusters)
        clusters = new_assignments[i] + 1;

    // define new clustering
    clusterings[way] = clustering_t(clusters, cluster_t());
    for(int i = 0; i < (int)new_assignments.size(); i++)
      clusterings[way][new_assignments[i]].push_back(i);

    return true;
  }

  bool Multiclustering::optimize(int way, int old_cluster, int index,
    vector<vector<counts_t> > & units_signatures,
    vector<vector<counts_t> > & clusters_signatures,
    vector<int> & new_assignments)
  // Library facilities used: assert
  {
    int & units = data->matrix.dimensions[way];
    assert(index < units);

    int unit = clusterings[way][old_cluster][index];
    int blocks = blocking_size(way);

    // search best cluster
    int new_cluster = -1;
    double new_cost = DBL_MAX;
    for(int cluster = 0; cluster != int(clusterings[way].size()); ++cluster)
    {
      // cost of encoding unit in way cluster
      double cluster_cost = 0;
      for(int block = 0; block != blocks; ++block)
        cluster_cost += hoffman_coding
        (units_signatures[unit][block], clusters_signatures[cluster][block]);

      // keep track of best cluster
      if(cluster_cost < new_cost)
      {
        new_cost = cluster_cost;
        new_cluster = cluster;
      }
    }

    // re-assign
    new_assignments[unit] = new_cluster;

    // feedback
    return (new_cluster != old_cluster);
  }
}
